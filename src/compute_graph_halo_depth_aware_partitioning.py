##
 #  Project: HAMPtOn
 #
 #  File: compute_graph_halo_depth_aware_partitioning.py
 #  Created: Jul 24, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys
import getopt
import argparse
import subprocess
import numpy as np
import random

from sets import Set
from operator import add, mul

from graph import Graph


MAJOR_VERSION=0
MINOR_VERSION=9


TEMPERATURE = 1e-3      ## for monte-carlo process convergence
A = 1.                  ## halo cells vs. local cells compute ratio
B = 1.                  ## communication vs. computation ratio: nhalocells / (ncellsvertical)
MAX_DEPTH = 40


class Partitioner:
  
  def __init__(self):
    self.ap = argparse.ArgumentParser()
    self.ap.add_argument('-i','--input', help='input graph file', required=True, type=file)
    self.ap.add_argument('-m','--metis', help='metis binary path', required=True)
    self.ap.add_argument('-w','--weight', help='input graph weight/depth file', required=False, type=file, default=None)
    self.ap.add_argument('-p','--partinit', help='initial partitioning file', required=False, type=file, default=None)
    self.ap.add_argument('-o','--output', help='output file prefix', required=True)
    self.ap.add_argument('-y','--nlayers', help='number of halo layers', required=False, default=0, choices=[0, 1, 2, 3], type=int)
    self.ap.add_argument('-n','--nparts', help='number of partitions', required=True, type=int)
    self.ap.add_argument('-t','--niters', help='number of iterations', required=False, default=10, type=int)
    self.ap.add_argument('-v', '--verbose', help='be verbose', required=False, default=0, action='count')
    self.ap.add_argument('-V', '--version', action='version', version='%(prog)s '+str(MAJOR_VERSION)+'.'+str(MINOR_VERSION))
    self.args = ap.parse_args()
    
    self.graph = Graph()


def compute_node_weights(nparts, parts, depths, halos, pneighbors, weights):
  ## compute num halo cells
  nhalocells = {}
  for p in range(0, nparts):
    nhalocells[p] = 0
  for p, hcells in halos.iteritems():
    nhalocells[p] = len(hcells)
  ## compute num cells
  ncells = {}
  for p in range(0, nparts):
    ncells[p] = 0
  for node, part in parts.iteritems():
    ncells[part] += 1
  ## compute num total cells
  ntotcells = map(add, ncells.values(), nhalocells.values())
  ## compute num neighbors
  nneighbors = {}
  for p in range(0, nparts):
    nneighbors[p] = 0
  for p, n in pneighbors.iteritems():
    nneighbors[p] = len(n)
  #print "*** nhalocells: ", nhalocells
  #print "*** ncells: ", ncells
  #print "*** total cells: ", ntotcells
  ## print "*** total cells stats: ",
  ## print_stats(ntotcells)
  ## compute cost due to computations
  local_comps = {}
  halo_comps = {}
  cell_comps = {}
  hcell_comps = {}
  compute_model_computation_cost(ncells, nhalocells, nparts, parts, depths, halos, local_comps, halo_comps, cell_comps, hcell_comps)
  tot_comps = {}
  for p in range(0, nparts):
    tot_comps[p] = local_comps[p] + halo_comps[p]
  ## compute cost due to communications
  tot_comms = {}
  compute_model_communication_cost(nhalocells, nneighbors, nparts, halos, depths, tot_comps, tot_comms)
  ## calculate the weights to add to each partition cell due to halo cells
  part_add_weights = {}
  for p in range(0, nparts):
    part_add_weights[p] = float(halo_comps[p] + tot_comms[p]) / ncells[p]
  for node, part in parts.iteritems():
    weights[node] = int((cell_comps[node] + part_add_weights[part]) * 100)


def compute_model_computation_cost(ncells, nhcells, nparts, parts, depths, halos, local_comps, halo_comps, cell_comps, hcell_comps):
  a = A
  f = F(nparts)
  for p in range(0, nparts):
    local_comps[p] = 0.
    halo_comps[p] = 0.
  for n, p in parts.iteritems():
    cell_comps[n] = depths[n] / f
    hcell_comps[n] = a * depths[n] / f
    local_comps[p] += cell_comps[n]
    if n in halos[p]: halo_comps[p] += hcell_comps[n]


def compute_model_communication_cost(nhcells, nneighbors, nparts, halos, depths, tot_comps, tot_comms):
  f = F(nparts)
  max_comp = max(tot_comps.values())
  ## all 40 layers are communicated
  for p in range(0, nparts):
    tot_comms[p] = B * ((40. * nhcells[p] / (nneighbors[p] * f)) + (max_comp - tot_comps[p]))


def calculate_parts_total_weights(nparts, parts, weights, pweights):
  for p in range(0, nparts):
    pweights[p] = 0
  for node in range(0, len(weights)):
    part = parts[node]
    pweights[part] += weights[node]


def print_stats(arr):
  pdarr = np.array(arr)
  asum = pdarr.sum()
  amin = pdarr.min()
  amax = pdarr.max()
  amean = pdarr.mean()
  astd = pdarr.std()
  print '\n++ Total weight: ' + str(asum)
  print '++ Min: ' + str(amin) + ', Max: ' + str(amax) + ', Mean: ' + str(amean) + ', Std: ' + str(astd)
  print '++ Imbalance: ' + str(1 - (float(amin) / amax))


def write_weighted_hypergraph(filename, data, weights, nnodes, nnets, npins):
  if nnodes != nnets: print 'warning: nnodes and nnets do not match: %d, %d\n' % (nnodes, nnets)
  ff = open(filename, 'w')
  header = '0 %d %d %d 1\n' % (nnodes, nnets, npins)
  ff.write(header)
  for net, pins in sorted(data.iteritems()):
    record = ''
    for p in pins: record += str(p) + ' '
    record += '\n'
    ff.write(record)
  record = ''
  for w in weights: record += str(w) + ' '
  record += '\n'
  ff.write(record)
  ff.close()


def write_parts(filename, parts):
  ff = open(filename, 'w')
  for node, part in sorted(parts.iteritems()):
    record = str(part) + '\n'
    ff.write(record)
  ff.close()



## the main part


##
## argument parsing
##

ap = argparse.ArgumentParser()
ap.add_argument('-i','--input', help='input graph file', required=True, type=file)
ap.add_argument('-m','--metis', help='metis binary path', required=True)
ap.add_argument('-w','--weight', help='input graph weight/depth file', required=False, type=file, default=None)
ap.add_argument('-p','--partinit', help='initial partitioning file', required=False, type=file, default=None)
ap.add_argument('-o','--output', help='output file prefix', required=True)
ap.add_argument('-y','--nlayers', help='number of halo layers', required=False, default=0, choices=[0, 1, 2, 3], type=int)
ap.add_argument('-n','--nparts', help='number of partitions', required=True, type=int)
ap.add_argument('-t','--niters', help='number of iterations', required=False, default=10, type=int)
ap.add_argument('-v', '--verbose', help='be verbose', required=False, default=0, action='count')
ap.add_argument('-V', '--version', action='version', version='%(prog)s '+str(MAJOR_VERSION)+'.'+str(MINOR_VERSION))

args = ap.parse_args()

graphfile = args.input
weightfile = args.weight
partfile = args.partinit
outfile = args.output
nparts = args.nparts
niters = args.niters
nlayers = args.nlayers
verbose = args.verbose


mypart = Partitioner()


sys.exit(2)

##
## the first run
##

data = {}  ## mapping { net -> list of nodes in net }
parts = {}  ## mapping { node -> partition number }
weights = {}  ## mapping { node -> cell cepth }

## the obtained data and graph always have base index as 0
nnodes, nnets, npins = read_partitioned_graph(graphfile, partfile, weightfile, data, weights, parts)
halos = {}
construct_parts_halos(data, nparts, parts, nlayers, halos)
pneighbors = {}
compute_part_neighbors(halos, nparts, parts, pneighbors)

## compute weights  

weights = map(int, np.ones(nnodes))
compute_node_weights(nparts, parts, depths, halos, pneighbors, weights)
#weights = map(weight_factor, weights)
totpweights = {}
calculate_parts_total_weights(nparts, parts, weights, totpweights)
## print "*** total part weights: ", totpweights
tmpoutfile = '.tmp.hypergraph'
write_weighted_hypergraph(tmpoutfile, data, weights, nnodes, nnets, npins)

infile = tmpoutfile
partfile = tmpoutfile + '.part.' + str(nparts)

err = 1e20  ## something big
random.seed()
logfile = open('.tmp.log', 'w')

## initial partitioning
final_parts = parts
final_weights = weights

## print "*** initial parts: ", parts
totpweights = {}
calculate_parts_total_weights(nparts, final_parts, weights, totpweights)
## print "*** initial total part weights: ", totpweights
print "** Initial partitioning weight stats: ",
print_stats(totpweights.values())

for i in range(0, niter):
  ## perform partitioning

  print "== Iteration", i, "=="
  seed = random.randint(1, 400000)
  cmd = patoh + ' ' + tmpoutfile + ' ' + str(nparts)
  #opts = 'UM=O PQ=S BO=C PA=3 RA=6 SD=' + str(seed)
  opts = 'SD=' + str(seed)
  cmd += ' ' + opts
  ## print cmd
  try:
    ret = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=logfile)
  except subprocess.CalledProcessError:
    print 'error: patoh failed!'
    sys.exit(2)
  else:
    if ret != 0:
      print 'error: something bad happened while running patoh!'
      sys.exit(2)
  parts = {}
  read_partitioning(partfile, parts)

  ## compute new halos based on new partitioning and neighbors
  halos = {}
  construct_parts_halos(data, nparts, parts, nlayers, halos)
  pneighbors = {}
  compute_part_neighbors(halos, nparts, parts, pneighbors)

  ## compute weights
  weights = map(int, np.ones(nnodes))
  compute_node_weights(nparts, parts, depths, halos, pneighbors, weights)
  #weights = map(weight_factor, weights)

  totpweights = {}
  calculate_parts_total_weights(nparts, parts, weights, totpweights)

  newerr = 1. - (float(min(totpweights.values())) / max(totpweights.values()))
  differr = err - newerr
  if differr > 0: accept = True
  else:
    prob = np.exp(differr / TEMPERATURE)
    ## print '*** Acceptance probability: ' + str(prob)
    if random.random() < prob: accept = True
    else: accept = False;
  
  ## print "*** ERROR VALUE: " + str(newerr) + " [ " + str(err) + ' ]'
  if not accept:
    ## print '*** rejecting partitioning'
    ## print "*** rejected part weights: ", totpweights
    ## print "*** rejected part weight stats: ",
    ## print_stats(totpweights.values())
    continue

  err = newerr
  ## print '*** accepting partitioning'  
  final_parts = parts
  final_weights = weights

  ## print "*** total part weights: ", totpweights
  ## print "*** part weight stats: ",
  ## print_stats(totpweights.values())

  if i == niter: break

  ## write new weighted hypergraph
  write_weighted_hypergraph(tmpoutfile, data, final_weights, nnodes, nnets, npins)

## write out the final partitioning. writing the hypergraph file is not needed ...
totpweights = {}
calculate_parts_total_weights(nparts, final_parts, final_weights, totpweights)
partoutfile = outfile + '.part.' + str(nparts)
write_parts(partoutfile, final_parts)

## print "*** final total part weights: ", totpweights
print "*** final part weight stats: ",
print_stats(totpweights.values())

logfile.close()
