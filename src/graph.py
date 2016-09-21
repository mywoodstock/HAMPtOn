##
 #  Project:
 #
 #  File: graph.py
 #  Created: Sep 18, 2016
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 #
 #  Metis graph file format:
 #    numnodes numedges format numconstraints
 #    s w1 w2 ... wncon v1 e1 v2 e2 ... vk ek     %% vertices adjacent to this vertex
 #    ...
 #  format = XYZ, Z = has edge weights, Y = has vertex weights, X = has vertex sizes
 #  numconstriaints = number of vertex weights (constraints)
 ##

import sys
import numpy as np

from sets import Set

class Graph:

  def __init__(self):
    self.nodelist = {}     ## map of "nodeid => { 'neighbors': nodelist, 'partition': partition }"
    self.partitions = {}   ## map of "partition => nodelist"
    self.total_num_nodes = 0
    self.total_num_edges = 0
  
  def num_nodes(self):
    return len(self.nodelist)

  def num_edges(self):
    return self.total_num_edges

  def neighbor_nodes(self, node):
    if node not in self.nodelist: return []
    return self.nodelist[node]

  def node_partition(self, node):
    if node not in self.nodelist: return -1
    return self.nodelist[node]['partition']
  
  def node_weight(self, node):
    if node not in self.nodelist: return 0
    return self.nodelist[node]['weight']
  
  def nodes_in_partition(self, part):
    if part not in self.partitions: return []
    return self.partitions[part]['nodes']
  
  def num_nodes_in_partition(self, part):
    if part not in self.partitions: return 0
    return len(self.partitions[part])  

  def num_partition_neighbors(self, part):
    if part not in self.partitions: return 0
    return len(self.partitions[part]['neighbors'])

  def num_partitions(self):
    return len(self.partitions)
  
  def set_node_weights(self, weights):
    for nodeid in range(1,self.num_nodes()+1):
      self.nodelist[nodeid]['weight'] = weights[nodeid]
    
  def read_graph(self, infile):
    print 'reading graph...'
    booleanify = lambda x: True if x == '1' else False
    gg = open(infile)
    line = gg.readline()            ## first line has num nodes, num edges, format, and num constraints
    words = line.strip().split()
    self.total_num_nodes = int(words[0])
    self.total_num_edges = int(words[1])
    has_vertex_sizes, has_vertex_weights, has_edge_weights = False, False, False
    num_vertex_weights = 0
    if len(words) > 2:
      digits = [ w for w in words[2] ]
      if len(digits) != 3: sys.exit(2)
      [ has_vertex_sizes, has_vertex_weights, has_edge_weights ] = map(booleanify, digits)
      if len(words) > 3: num_vertex_weights = int(words[3])
      else: num_vertex_weights = 1
    nodeid = 0
    while True:
      line = gg.readline()
      if not line: break
      nodeid += 1
      words = map(int, line.strip().split())
      if nodeid not in self.nodelist: self.nodelist[nodeid] = {}
      index = 0
      if has_vertex_sizes:
        self.nodelist[nodeid]['size'] = words[index]
        index += 1
      if has_vertex_weights:
        self.nodelist[nodeid]['weight'] = words[index:num_vertex_weights+index]
        index += num_vertex_weights
      if has_edge_weights:
        neighbors = words[index::2]
        edgeweights = words[index+1::2]
      else:
        neighbors = words[index:]
        edgeweights = [ 1 ] * len(neighbors)
      self.nodelist[nodeid]['neighbors'] = zip(neighbors, edgeweights)
      if 'partition' not in self.nodelist[nodeid]: self.nodelist[nodeid]['partition'] = 0   ## initialize
    gg.close()
    if nodeid != self.total_num_nodes:
			print 'error: number of nodes mismatch in graph file'
			sys.exit(2)
  
  def read_weights(self, wfile=None, weights=None, sfile=None, sizes=None):
    print 'initializing weights...'
    ## node weights: local weight, total (local+halo) weight (computation)
    ## node size (communication volume)
    if wfile is not None:
      ff = open(wfile)
      nodeid = 0
      while True:
        line = ff.readline()
        if not line: break
        words = line.strip().split()
        if len(words) > 0:
          for word in words:
            nodeid += 1
            self.nodelist[nodeid]['weight'] = [ float(word), 0. ]
      ff.close()
    else:
      if weights is None: weights = 1.
      for nodeid in range(1, self.num_nodes()+1):
        self.nodelist[nodeid]['weight'] = [ weights, 0. ]
    if sfile is not None:
      print 'uh-oh: reading sizes from a file has not been implemented!'
    else:
      if sizes is None: sizes = 1.
      for nodeid in range(1, self.num_nodes()+1):
        self.nodelist[nodeid]['size'] = sizes

  def read_partitions(self, infile):
    print 'reading partitions...'
    ff = open(infile)
    nodeid = 0
    while True:
      line = ff.readline()
      if not line: break
      parts = map(int, line.strip().split())
      for part in parts:
        nodeid += 1
        if part not in self.partitions:
          self.partitions[part] = {}
          self.partitions[part]['nodes'] = set()
          self.partitions[part]['halo_nodes'] = set()
          self.partitions[part]['neighbors'] = set()
        self.partitions[part]['nodes'].add(nodeid)
        if nodeid not in self.nodelist: self.nodelist[nodeid] = {}
        self.nodelist[nodeid]['partition'] = part
    ff.close()

  def get_halo_candidates(self, candidates, node, nlayers):
		for (neighbor, eweight) in self.nodelist[node]['neighbors']:
			candidates.add(neighbor)
			if nlayers > 0: self.get_halo_candidates(candidates, neighbor, nlayers-1)

  def compute_halos(self, nlayers):
    print "computing halos ..."
    if nlayers < 1: return
    for part in range(0, self.num_partitions()):
			candidates = set()
			for node in self.nodes_in_partition(part):
				self.get_halo_candidates(candidates, node, nlayers)
			self.partitions[part]['halo_nodes'] = list(candidates - set(self.partitions[part]['nodes']))
	
  def compute_part_neighbors(self):
		for part, pdata in self.partitions.iteritems():
			for hnode in pdata['halo_nodes']:
				self.partitions[part]['neighbors'].add(self.nodelist[hnode]['partition'])

  def compute_local_weights(self):
    local_weights = {}
    for p in range(0, self.num_partitions()):
      pnodes = self.partitions[p]['nodes']
      local_weights[p] = np.array([self.nodelist[nodeid]['weight'][0] for nodeid in pnodes]).sum()
    return local_weights

  def compute_halo_weights(self):
    halo_weights = {}
    for p in range(0, self.num_partitions()):
      pnodes = self.partitions[p]['halo_nodes']
      halo_weights[p] = np.array([self.nodelist[nodeid]['weight'][0] for nodeid in pnodes]).sum()
    return halo_weights

  def update_local_weights(self, weights):
    for p in range(0, self.num_partitions()):
      pnodes = self.partitions[p]['nodes']
      hweight = float(weights[p]) / len(pnodes)
      for nodeid in pnodes:
        self.nodelist[nodeid]['weight'][1] = hweight + self.nodelist[nodeid]['weight'][0]

  def compute_halo_volumes(self):
    halo_volumes = {}
    for p in range(0, self.num_partitions()):
      pnodes = self.partitions[p]['halo_nodes']
      sizes = np.array([self.nodelist[nodeid]['size'] for nodeid in pnodes])
      halo_volumes[p] = sizes.sum()
    return halo_volumes

  def print_statistics(self):
    num_parts = self.num_partitions()
    temp_part_num_cells = []
    temp_part_num_halo_cells = []
    temp_part_num_neighbors = []
    for part, partdata in self.partitions.iteritems():
      temp_part_num_cells.append(len(partdata['nodes']))
      temp_part_num_halo_cells.append(len(partdata['halo_nodes']))
      temp_part_num_neighbors.append(len(partdata['neighbors']))
    part_num_nodes = np.array(temp_part_num_cells)
    part_num_halos = np.array(temp_part_num_halo_cells)
    part_num_total = part_num_nodes + part_num_halos
    part_num_neighbors = np.array(temp_part_num_neighbors)
    
    part_weight_nodes = {}
    part_weight_halos = {}
    part_weight_total = {}
    part_volume_halos = {}
    for part, pdata in self.partitions.iteritems():
      part_weight_nodes[part] = 0
      for nid in pdata['nodes']:
        part_weight_nodes[part] += self.nodelist[nid]['weight'][0]
      part_weight_halos[part] = 0
      part_volume_halos[part] = 0
      for nid in pdata['halo_nodes']:
        part_weight_halos[part] += self.nodelist[nid]['weight'][0]
        part_volume_halos[part] += self.nodelist[nid]['size'] / len(pdata['neighbors'])
      part_weight_total[part] = part_weight_nodes[part] + part_weight_halos[part]
    part_weight_nodes = np.array(part_weight_nodes.values())
    part_weight_halos = np.array(part_weight_halos.values())
    part_weight_total = np.array(part_weight_total.values())
    part_volume_halos = np.array(part_volume_halos.values())
    
    print ("**        num nodes:"
				+ " min = " + str(part_num_nodes.min())
				+ " max = " + str(part_num_nodes.max())
				+ " mean = " + str(part_num_nodes.mean())
				+ " std-dev = "	+ str(part_num_nodes.std())
				+ " imbalance = " + str(1.0 - float(part_num_nodes.min())/part_num_nodes.max()))
    print ("**   num halo nodes:"
				+ " min = " + str(part_num_halos.min())
				+ " max = " + str(part_num_halos.max())
				+ " mean = " + str(part_num_halos.mean())
				+ " std-dev = " + str(part_num_halos.std())
				+ " imbalance = " + str(1.0 - float(part_num_halos.min())/part_num_halos.max()))
    print ("**  total num nodes:"
				+ " min = " + str(part_num_total.min())
				+ " max = " + str(part_num_total.max())
				+ " mean = " + str(part_num_total.mean())
				+ " std-dev = " + str(part_num_total.std())
				+ " imbalance = " + str(1.0 - float(part_num_total.min())/part_num_total.max()))
    print ("**        neighbors:"
				+ " min = " + str(part_num_neighbors.min())
				+ " max = " + str(part_num_neighbors.max())
				+ " mean = " + str(part_num_neighbors.mean())
				+ " std-dev = " + str(part_num_neighbors.std())
				+ " imbalance = " + str(1.0 - float(part_num_neighbors.min())/part_num_neighbors.max()))
    
    print ("**     weight nodes:"
				+ " min = " + str(part_weight_nodes.min())
				+ " max = " + str(part_weight_nodes.max())
				+ " mean = " + str(part_weight_nodes.mean())
				+ " std-dev = "	+ str(part_weight_nodes.std())
				+ " imbalance = " + str(1.0 - float(part_weight_nodes.min())/part_weight_nodes.max()))
    print ("**     weight halos:"
				+ " min = " + str(part_weight_halos.min())
				+ " max = " + str(part_weight_halos.max())
				+ " mean = " + str(part_weight_halos.mean())
				+ " std-dev = "	+ str(part_weight_halos.std())
				+ " imbalance = " + str(1.0 - float(part_weight_halos.min())/part_weight_halos.max()))
    print ("**     weight total:"
				+ " min = " + str(part_weight_total.min())
				+ " max = " + str(part_weight_total.max())
				+ " mean = " + str(part_weight_total.mean())
				+ " std-dev = "	+ str(part_weight_total.std())
				+ " imbalance = " + str(1.0 - float(part_weight_total.min())/part_weight_total.max()))

    print ("**     volume halos:"
				+ " min = " + str(part_volume_halos.min())
				+ " max = " + str(part_volume_halos.max())
				+ " mean = " + str(part_volume_halos.mean())
				+ " std-dev = "	+ str(part_volume_halos.std())
				+ " imbalance = " + str(1.0 - float(part_volume_halos.min())/part_volume_halos.max()))
    

  def print_detailed_statistics(self):
		print "** partitioning sizes:"
		for part, pdata in self.partitions.iteritems():
			nnodes = len(pdata['nodes'])
			nhnodes = len(pdata['halo_nodes'])
			nneighbors = len(pdata['neighbors'])
			print ("    partition " + str(part) + " :: nodes = " + str(nnodes)
					+ " :: halo nodes = " + str(nhnodes)
					+ " :: total nodess = " + str(nhnodes + nnodes)
					+ " :: neighbors = " + str(nneighbors))
		self.print_statistics()

  def printall(self):
    print self.nodelist
    print self.partitions

  def save_graph(self, filename):
    flatten_neighbors = lambda (v, e): str(v) + '\t' + str(int(e))
    ff = open(filename, 'w')
    header = '%d %d 111 1\n' % (self.num_nodes(), self.num_edges())
    ff.write(header)
    for nodeid in range(1,self.num_nodes()+1):
      rec = ''
      rec += str(int(self.nodelist[nodeid]['size'])) + '\t'
      rec += str(int(self.nodelist[nodeid]['weight'][1])) + '\t'                    ## save only the second weight
      rec += '\t'.join(map(flatten_neighbors, self.nodelist[nodeid]['neighbors']))
      rec += '\n'
      ff.write(rec)
    ff.close()

  def save_partitioning(self, filename):
    ff = open(filename, 'w')
    part_list = [ str(self.nodelist[n]['partition']) for n in range(1, self.num_nodes()+1) ]
    ff.write('\n'.join(part_list) + '\n')
    ff.close()

  def save_partition_info(self, oprefix):
		## save halos and partition neighbor graph
		hfilename = oprefix + '_partition_halos.info'
		ff = open(hfilename, 'w')
		p = str(self.num_partitions()) + '\n'
		ff.write(p)
		for part, pdata in sorted(self.partitions.iteritems()):
			phalo = '\t'.join(map(str, pdata['halo_nodes'])) + '\n'
			ff.write(phalo)
		ff.close()
		pfilename = oprefix + '_partition_graph.info'
		ff = open(pfilename, 'w')
		p = str(self.num_partitions()) + '\n'
		ff.write(p)
		for part, pdata in sorted(self.partitions.iteritems()):
			pneighbors = '\t'.join(map(str, pdata['neighbors'])) + '\n'
			ff.write(pneighbors)
		ff.close()
    
  def read_partition_info(self, oprefix):
		hfilename = oprefix + '_partition_halos.info'
		print 'reading halos from', hfilename, '...'
		ff = open(hfilename, 'r')
		nparts = int(ff.readline().strip().split()[0])
		part = 0
		while True:
			line = ff.readline()
			if not line: break
			self.partitions[part]['halo_nodes'] = map(int, line.strip().split())
			part += 1
		ff.close()
		if nparts != part:
			print 'error: mismatching number of partitions'
			sys.exit(2)
		pfilename = oprefix + '_partition_graph.info'
		print 'reading partition graph from', pfilename, '...'
		ff = open(pfilename, 'r')
		nparts = int(ff.readline().strip().split()[0])
		part = 0
		while True:
			line = ff.readline()
			if not line: break
			self.partitions[part]['neighbors'] = map(int, line.strip().split())
			part += 1
		ff.close()
		if nparts != part:
			print 'error: mismatching number of partitions'
			sys.exit(2)