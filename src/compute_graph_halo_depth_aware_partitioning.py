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
    ##
    ## argument parsing
    ##
    ap = argparse.ArgumentParser()
    ap.add_argument('-i','--input', help='input graph file', required=True)
    ap.add_argument('-m','--metis', help='metis binary path', required=True)
    ap.add_argument('-w','--weight', help='input graph weight/depth file', required=False, default=None)
    ap.add_argument('-p','--partinit', help='initial partitioning file', required=False, default=None)
    ap.add_argument('-o','--output', help='output file prefix', required=True)
    ap.add_argument('-y','--nlayers', help='number of halo layers', required=False, default=0, choices=[0, 1, 2, 3], type=int)
    ap.add_argument('-n','--nparts', help='number of partitions', required=True, type=int)
    ap.add_argument('-t','--niters', help='number of iterations', required=False, default=10, type=int)
    ap.add_argument('-v', '--verbose', help='be verbose', required=False, default=0, action='count')
    ap.add_argument('-V', '--version', action='version', version='%(prog)s '+str(MAJOR_VERSION)+'.'+str(MINOR_VERSION))
    args = ap.parse_args()
    self.inputconfig = {
      'metis': args.metis,
      'graphfile': args.input,
      'weightfile': args.weight,
      'partfile': args.partinit,
      'outprefix': args.output,
      'nparts': args.nparts,
      'niters': args.niters,
      'nlayers': args.nlayers,
      'verbose': args.verbose
    }
    self.logfile = open('hampton.log', 'a')
    ##
    ## initialize graph
    ##
    self.graph = Graph()    
    self.graph.read_graph(self.inputconfig['graphfile'])
    if self.inputconfig['weightfile'] is not None:
      self.graph.read_weights(wfile=self.inputconfig['weightfile'],sizes=MAX_DEPTH)
    init_partfile = self.perform_init_partitioning()
    self.graph.read_partitions(init_partfile)
    ##
    ## compute halos and neighbors
    ##
    self.graph.compute_halos(self.inputconfig['nlayers'])
    self.graph.compute_part_neighbors()
    ##
    ## misc stuff
    ##
    self.num_cells = self.graph.num_nodes()
    self.num_partitions = self.inputconfig['nparts']
    if self.num_partitions != self.graph.num_partitions():
      print 'error: number of partitions is not right'
      sys.exit(2)
    ##
    ## display graph details
    ##
    self.graph.printall()
    self.graph.print_detailed_statistics()
  
  def perform_init_partitioning(self):
    if self.inputconfig['partfile'] is not None:
      print 'info: using given initial partitioning'
      return self.inputconfig['partfile']
    print 'info: performing initial partitioning'
    seed = random.randint(1, 400000)
    cmd = self.inputconfig['metis'] + ' ' + self.inputconfig['graphfile'] + ' ' + str(self.inputconfig['nparts'])
    opts = ''
    opts = '-seed=' + str(seed)
    cmd += ' ' + opts
    print 'executing', cmd
    try:
      ret = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=self.logfile)
    except subprocess.CalledProcessError:
      print 'error: metis failed! check log for details.'
      sys.exit(2)
    else:
      if ret != 0:
        print 'error: something bad happened while running metis! check the log for details'
        sys.exit(2)
    init_partfile = self.inputconfig['graphfile'] + '.part.' + str(self.inputconfig['nparts'])
    return init_partfile

  def compute_and_update_weights(self):
    ## compute 'local' computation weights for each partition
    local_weights = self.graph.compute_local_weights()
    halo_weights = self.graph.compute_halo_weights()
    total_weights = { p: local_weights[p] + halo_weights[p] for p in range(0, self.graph.num_partitions())}
    self.graph.update_local_weights(halo_weights)
    ## compute 'halo' computation weight for each partition
    halo_volumes = self.graph.compute_halo_volumes()
    ## compute 'halo' communication volume (total for each partition)
    ## compute 'halo' communication volume (total send volume and per neighbor, total receive volume and per neighbor)
    print 'local weights:', local_weights
    print 'halo weights :', halo_weights
    print 'total weights:', total_weights
    print 'halo volumes :', halo_volumes

  def write_partitioning(self):
    filename = self.inputconfig['outprefix'] + '.info.part.' + str(self.inputconfig['nparts'])
    self.graph.save_partitioning(filename)

  def write_weighted_graph(self):
    filename = 'tmp.' + self.inputconfig['outprefix'] + '.info'
    self.graph.save_graph(filename)


## the main part

hampton = Partitioner()
hampton.compute_and_update_weights()
hampton.write_weighted_graph()
hampton.write_partitioning()

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
