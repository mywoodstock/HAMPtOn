##
 #  Project:
 #
 #  File: analyzepartitions.py
 #  Created: Oct 10, 2016
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt
import subprocess
import numpy as np
import pandas as pd
import random

from sets import Set
from operator import add, mul


def parse_arguments(argv):
  try:
    opts, args = getopt.getopt(argv, 'g:p:d:o:y:n:')
    if len(opts) != 6:
      raise getopt.GetoptError('Give arguments')
  except getopt.GetoptError:
    print 'analyzepartitions.py -g <graphfile> -p <partitionfile> -d <depthfile> -o <outputfile> -y <nhalolayers> -n <nparts>'
    sys.exit(2)
  for opt, arg in opts:
    if opt in ('-g'): graphfile = arg
    elif opt in ('-o'): outputfile = arg
    elif opt in ('-p'): partfile = arg
    elif opt in ('-d'): depthfile = arg
    elif opt in ('-n'): nparts = int(arg)
    elif opt in ('-y'): nhalos = int(arg)
  return graphfile, partfile, depthfile, outputfile, nparts, nhalos


def read_partitioning(partfile, nodes):
  count = 0
  ff = open(partfile)
  while True:
    line = ff.readline()
    if not line: break
    words = map(int, line.strip().split())
    for word in words:
      count += 1
      nodes[count] = word
  ff.close()


def read_depths(depthfile, depths):
  node = 0
  ff = open(depthfile)
  while True:
    line = ff.readline()
    if not line: break
    words = map(int, line.strip().split())
    for word in words:
      node += 1
      depths[node] = word
  ff.close()


def group_partitions(nodes, parts):
  for node, part in nodes.iteritems():
    if part not in parts: parts[part] = []
    parts[part].append(node)


def read_partitioned_graph(graphfile, partfile, depthfile, graph, nodes, depths, parts):
  ## first read the graph
  ff = open(graphfile)
  line = ff.readline()    ## first line has [ #nodes, #edges ]
  if not line: return
  words = map(int, line.strip().split())
  if len(words) != 2:
    print 'error: invalid header in graph file'
    sys.exit(2)
  nnodes = words[0]
  nedges = words[1]
  node = 0
  while True:
    line = ff.readline()
    if not line: break
    node += 1
    neighbors = map(int, line.strip().split())
    graph[node] = neighbors
  ff.close()
  if node != nnodes:
    print 'error: mismatch in number of nodes and lines in graph file'
    sys.exit(2)
  ## read the partition file
  read_partitioning(partfile, nodes)
  ## group nodes into partitions
  group_partitions(nodes, parts)
  ## read the depth file
  read_depths(depthfile, depths)
  return nnodes, nedges


def get_halo_candidates(candidates, visited, graph, node, nlayers):
  if nlayers < 1: return
  for neighbor in graph[node]:
    candidates.add(neighbor)
    if neighbor not in visited:
      visited.add(neighbor)
      get_halo_candidates(candidates, visited, graph, neighbor, nlayers-1)


def compute_halos(graph, nparts, parts, nlayers, halos):
  print "computing halos ..."
  if nlayers < 1: return
  for part in range(0, nparts):
    candidates = set()
    visited = set()
    for node in parts[part]:
      visited.add(node)
      get_halo_candidates(candidates, visited, graph, node, nlayers)
    halos[part] = list(candidates - set(parts[part]))


def compute_part_weights(nnodes, nparts, parts, halos, depths, outfile):
  part_nnodes = {}
  part_nhalos = {}
  part_nvertnodes = {}
  part_nverthalos = {}
  for p in range(0, nparts):
    part_nnodes[p] = len(parts[p])
    part_nhalos[p] = len(halos[p])
    nvertnodes = 0
    nverthalos = 0
    for n in parts[p]:
      nvertnodes += depths[n]
    for n in halos[p]:
      nverthalos += depths[n]
    part_nvertnodes[p] = nvertnodes
    part_nverthalos[p] = nverthalos
  return part_nnodes, part_nhalos, part_nvertnodes, part_nverthalos


def write_stats(outfile, nnodes, nparts, part_nnodes, part_nhalos, part_nvertnodes, part_nverthalos):
  ff = open(outfile, 'w')
  header = 'partition  ncells  nhcells  ntotcells  nvertcells nverthcells nverttotcells\n'
  ff.write(header)
  for p in range(0, nparts):
    record = '%d  %d  %d  %d  %d  %d  %d\n' % (p, part_nnodes[p], part_nhalos[p], part_nnodes[p]+part_nhalos[p], part_nvertnodes[p], part_nverthalos[p], part_nvertnodes[p]+part_nverthalos[p])
    ff.write(record)
  ff.close()


## the main part
## for a given partitioning, output:
##  { partition ncells nhcells ntotcells nvertcells nverthcells nverttotcells }

graphfile, partfile, depthfile, outfile, nparts, nlayers = parse_arguments(sys.argv[1:])

graph = {}  ## mapping { node -> list of neighbors of node }
nodes = {}  ## mapping { node -> partition number }
depths = {} ## mapping { node -> depth }
parts = {}  ## mapping { part -> list of nodes in part }
halos = {}  ## mapping { part -> list of nodes in part's halos }
nnodes, nedges = read_partitioned_graph(graphfile, partfile, depthfile, graph, nodes, depths, parts)
compute_halos(graph, nparts, parts, nlayers, halos)
part_nnodes, part_nhalos, part_nvertnodes, part_nverthalos = compute_part_weights(nnodes, nparts, parts, halos, depths, outfile)
write_stats(outfile, nnodes, nparts, part_nnodes, part_nhalos, part_nvertnodes, part_nverthalos)