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

class Graph:

  nodelist = {}             ## map of "nodeid => { 'neighbors': nodelist, 'partition': partition }"
  partitions = {}           ## map of "partition => nodelist"

  def __init__(self):
    self.partitions[0] = []
  
  def num_nodes(self):
    return len(self.nodelist)

  def neighbor_nodes(self, node):
    if node not in self.nodelist: return []
    return self.nodelist[node]

  def node_partition(self, node):
    if node not in self.nodelist: return -1
    return nodelist[node]['partition']
  
  def nodes_in_partition(self, part):
    if part not in self.partitions: return []
    return self.partitions[part]
  
  def num_nodes_in_partition(self, part):
    if part not in self.partitions: return 0
    return len(self.partitions[part])

  def num_partitions(self):
    return len(self.partitions)
    
  def read_graph(self, infile):
    print 'reading graph file', infile, '...'
    booleanify = lambda x: True if x == '1' else False
    gg = open(infile)
    line = gg.readline()            ## first line has num nodes, num edges, format, and num constraints
    words = line.strip().split()
    num_nodes = int(words[0])
    num_edges = int(words[1])
    has_vertex_sizes, has_vertex_weights, has_edge_weights = False, False, False
    num_vertex_weights = 0
    if len(words) > 2:
      digits = words[2].split('')
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
        self.nodelist[nodeid]['weights'] = words[index:num_vertex_weights+index]
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
    if nodeid != num_nodes:
			print 'error: number of nodes mismatch in graph file'
			sys.exit(2)
  
  def read_weights(self, infile):
    print 'reading weights file', infile, '...'
    ff = open(infile)
    nodeid = 0
    while True:
      line = ff.readline()
      if not line: break
      words = line.strip().split()
      if len(words) > 0:
        for word in words:
          nodeid += 1
          nodelist[nodeid]['weight'] = float(word)
    ff.close()

  def read_partitions(self, infile):
    print 'reading partitions file', infile, '...'
    pp = open(infile)
    nodeid = 0
    while True:
      line = pp.readline()
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
    pp.close()

  def get_halo_candidates(self, candidates, node, nlayers):
		for neighbor in self.nodelist[node]:
			candidates.add(neighbor)
			if nlayers > 0: self.get_halo_candidates(candidates, neighbor, nlayers-1)

  def compute_halos(self, nlayers):
		print "computing halos ..."
		for part in range(0, self.num_partitions()):
			candidates = set()
			for node in self.nodes_in_partition(part):
				self.get_halo_candidates(candidates, node, nlayers)
			self.partitions[part]['halo_nodes'] = list(candidates - set(self.partitions[part]['nodes']))
	
  def compute_part_neighbors(self):
		for part, pdata in self.partitions.iteritems():
			for hnode in pdata['halo_nodes']:
				pdata['neighbors'].add(nodelist[hnode]['partition'])

  def print_statistics(self):
		num_parts = self.num_partitions()
		temp_part_num_cells = []
		temp_part_num_halo_cells = []
		temp_part_num_neighbors = []
		for part, partdata in self.graph['partitions'].iteritems():
			temp_part_num_cells.append(len(partdata['cells']))
			temp_part_num_halo_cells.append(len(partdata['halo_cells']))
			temp_part_num_neighbors.append(len(partdata['neighbors']))
		part_num_cells = np.array(temp_part_num_cells)
		part_num_halo_cells = np.array(temp_part_num_halo_cells)
		part_num_neighbors = np.array(temp_part_num_neighbors)
		part_total_num_cells = part_num_cells + part_num_halo_cells
		print ("**            num cells:"
				+ " min = " + str(part_num_cells.min())
				+ " max = " + str(part_num_cells.max())
				+ " mean = " + str(part_num_cells.mean())
				+ " std-dev = "	+ str(part_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_cells.min())/part_num_cells.max()))
		print ("**       num halo cells:"
				+ " min = " + str(part_num_halo_cells.min())
				+ " max = " + str(part_num_halo_cells.max())
				+ " mean = " + str(part_num_halo_cells.mean())
				+ " std-dev = " + str(part_num_halo_cells.std())
				+ " imbalance = " + str(1.0 - float(part_num_halo_cells.min())/part_num_halo_cells.max()))
		print ("**      total num cells:"
				+ " min = " + str(part_total_num_cells.min())
				+ " max = " + str(part_total_num_cells.max())
				+ " mean = " + str(part_total_num_cells.mean())
				+ " std-dev = " + str(part_total_num_cells.std())
				+ " imbalance = " + str(1.0 - float(part_total_num_cells.min())/part_total_num_cells.max()))
		print ("**            neighbors:"
				+ " min = " + str(part_num_neighbors.min())
				+ " max = " + str(part_num_neighbors.max())
				+ " mean = " + str(part_num_neighbors.mean())
				+ " std-dev = " + str(part_num_neighbors.std())
				+ " imbalance = " + str(1.0 - float(part_num_neighbors.min())/part_num_neighbors.max()))
		nparts = len(self.graph['partitions'])

  def print_detailed_statistics(self):
		print "** partitions data:"
		for part, pdata in self.graph['partitions'].iteritems():
			ncells = len(pdata['cells'])
			nhcells = len(pdata['halo_cells'])
			nneighbors = len(pdata['neighbors'])
			print ("    partition " + str(part) + " :: cells = " + str(ncells)
					+ " :: halo cells = " + str(nhcells)
					+ " :: total cells = " + str(nhcells + ncells)
					+ " :: neighbors = " + str(nneighbors))
		self.print_statistics()

  def printall(self):
    print self.nodelist
    print self.partitions
  
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


if __name__ == '__main__':
  sys.exit(0)