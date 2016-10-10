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


TEMPERATURE = 1e-2        ## for monte-carlo process convergence
COMP_COMM_RATIO = 0.0     ## communication vs. computation ratio
COMP_SIZE_RATIO = 10.0    ## compute weight + COMP_SIZE_RATIO * compute size


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
    ap.add_argument('-d','--maxdepth', help='maximum depth of cells', required=False, default=1, type=int)
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
      'maxdepth': args.maxdepth,
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
      self.graph.read_weights(wfile=self.inputconfig['weightfile'],sizes=self.inputconfig['maxdepth'])
    self.store_node_depths()
    self.store_node_volumes()
    init_partfile = self.construct_init_partitioning()
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
    self.print_statistics()
  
  def construct_init_partitioning(self):
    if self.inputconfig['partfile'] is not None:
      print 'info: using given initial partitioning'
      return self.inputconfig['partfile']
    print 'info: performing initial partitioning'
    return self.construct_metis_partitioning(self.inputconfig['graphfile'], self.inputconfig['nparts'])

  def construct_metis_partitioning(self, graphfile, nparts):
    seed = random.randint(1, 400000)
    cmd = self.inputconfig['metis'] + ' ' + graphfile + ' ' + str(nparts)
    opts = '-seed=' + str(seed)
    cmd += ' ' + opts
    print 'executing', cmd, '...'
    try:
      ret = subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=self.logfile)
    except subprocess.CalledProcessError:
      print 'error: metis failed! check log for details.'
      sys.exit(2)
    else:
      if ret != 0:
        print 'error: something bad happened while running metis! check the log for details'
        sys.exit(2)
    partfile = graphfile + '.part.' + str(nparts)
    return partfile

  def store_node_depths(self):
    self.depths = { nodeid: self.graph.node_weight(nodeid) for nodeid in range(1,self.graph.num_nodes()+1) }

  def store_node_volumes(self):
    self.volumes = { nodeid: self.graph.node_size(nodeid) for nodeid in range(1,self.graph.num_nodes()+1) }

  def compute_and_update_weights(self, graph=None):
    if graph is None: graph = self.graph
    ## compute 'local' computation weights for each partition
    ## compute 'halo' computation weight for each partition
    local_weights = graph.compute_local_weights()
    halo_weights = graph.compute_halo_weights()
    total_weights = { p: local_weights[p] + halo_weights[p] for p in range(0, graph.num_partitions())}
    local_sizes = graph.compute_local_sizes()
    halo_sizes = graph.compute_halo_sizes()
    total_sizes = { p: local_sizes[p] + halo_sizes[p] for p in range(0, graph.num_partitions())}
    halo_total_weights = { p: halo_weights[p] + COMP_SIZE_RATIO * halo_sizes[p] for p in range(0, graph.num_partitions())}
    graph.update_local_weights(halo_total_weights, halo_mfactor=10, local_mfactor=10)
    ## compute 'halo' communication volume (total for each partition)
    ## compute 'halo' communication volume (total send volume and per neighbor, total receive volume and per neighbor)
    halo_volumes = graph.compute_halo_volumes()
    avg_halo_volumes = { p: halo_volumes[p] / graph.num_partition_neighbors(p) for p in range(0, graph.num_partitions()) }
    graph.update_local_volumes(halo_volumes)
    print 'local weights:', local_weights
    print 'halo weights :', halo_weights
    print 'total weights:', total_weights
    print 'halo volumes :', halo_volumes
    return total_weights, total_sizes, halo_volumes, avg_halo_volumes

  def construct_halo_aware_partitioning(self):
    err = 1e20    ## something big
    random.seed()
    temp_graphfile = 'tmp.' + self.inputconfig['outprefix'] + '.info'
    curr_graph = self.graph
    nparts = self.inputconfig['nparts']
    for i in range(0, self.inputconfig['niters']):
      print '## Running iteration', i, '...'
      test_partfile = self.construct_metis_partitioning(temp_graphfile, self.inputconfig['nparts'])
      test_graph = Graph()    
      test_graph.read_graph(temp_graphfile)
      test_graph.read_partitions(test_partfile)
      test_graph.set_node_weights(self.depths)
      test_graph.set_node_sizes(self.volumes)
      test_graph.compute_halos(self.inputconfig['nlayers'])
      test_graph.compute_part_neighbors()
      compute_weight_cost, compute_size_cost, halo_cost, avg_halo_cost = self.compute_and_update_weights(test_graph)
      compute_cost = { p: compute_weight_cost[p] + COMP_SIZE_RATIO * compute_size_cost[p] for p in range(0, nparts) }
      total_cost = { p: compute_cost[p] + COMP_COMM_RATIO * halo_cost[p] for p in range(0, nparts) }
      print 'compute cost:', compute_cost
      print 'halo cost:', halo_cost
      print 'total cost:', total_cost
      newerr = 1. - (float(min(total_cost.values())) / max(total_cost.values()))
      print '** Current error value:', newerr
      print '** Previous error value:', err
      differr = err - newerr
      if differr > 0: accept = True
      else:
        prob = np.exp(differr / TEMPERATURE)
        print '** Acceptance probability: ' + str(prob)
        if random.random() < prob: accept = True
        else: accept = False;
      if not accept:
        print '## Partitioning rejected'
        continue
      print '## Partitioning accepted'
      err = newerr
      curr_graph = test_graph
      self.write_weighted_graph(temp_graphfile, curr_graph)
    self.write_partitioning(graph=curr_graph)
    self.print_statistics(graph=curr_graph)

  def write_partitioning(self, filename=None, graph=None):
    if filename is None: filename = self.inputconfig['outprefix'] + '.info.part.' + str(self.inputconfig['nparts'])
    if graph is None: graph = self.graph
    graph.save_partitioning(filename)

  def write_weighted_graph(self, filename=None, graph=None):
    if filename is None: filename = 'tmp.' + self.inputconfig['outprefix'] + '.info'
    if graph is None: graph = self.graph
    graph.save_graph(filename)

  def print_statistics(self, graph=None):
    if graph is None: graph = self.graph
    # graph.printall()
    graph.print_detailed_statistics()


##
## the main part
##

if __name__ == '__main__':
  hampton = Partitioner()
  hampton.compute_and_update_weights()
  hampton.write_weighted_graph()
  hampton.write_partitioning()
  hampton.print_statistics()
  
  hampton.construct_halo_aware_partitioning()  
