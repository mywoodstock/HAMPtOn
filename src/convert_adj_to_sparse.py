##
 #  Project:
 #
 #  File: convert_adj_to_sparse.py
 #  Created: Mar 08, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 #
 #	Description: Convert graph.ino style adjacency graph data into sparse matrix representation
 #					such as tuples: (row, col, val), and add the MM header.
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print "convert_adj_to_sparse.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
	return inputfile, outputfile


def read_adj_file(filename, data):
	ff = open(filename)
	head = ff.readline()		# this has num nodes and num edges
	while head.startswith("%"): head = ff.readline()
	words = map(int, head.strip().split())
	if len(words) == 2:		## nodes edges
		nodes = words[0]
		edges = words[1]
		index = 1
	elif len(words) == 3:	## nx ny nnz
		nodes = words[0]
		nodes2 = words[1]
		edges = words[2]
		index = 1
		if nodes != nodes2:
			print "error: either the matrix is not square, or the header is malformed"
			sys.exit(2)
	elif len(words) == 4:	## patoh: index nx ny nnz
		nodes = words[1]
		nodes2 = words[2]
		edges = words[3]
		index = words[0]
		if nodes != nodes2:
			print "error: either the matrix is not square, or the header is malformed"
			sys.exit(2)
	else:
		print "error: unknown header type!"
		sys.exit(2)
	count = index
	offset = 1 - index		## assuming index is either 0 or 1
	nedges = 0
	while True:
		line = ff.readline()
		if not line: break
		data[count + offset] = []
		for word in line.strip().split():
			data[count + offset] += [int(word) + offset]
			nedges += 1
		count += 1
	ff.close()
	if count != nodes + index:
		print "error: mismatching number of nodes and data in the file: %d, %d" % (nodes, count)
		sys.exit(2)
	if nedges != edges:
		print "warning: mismatching number of edges and data in the file: %d, %d" % (edges, nedges)
	return nodes, nedges


def convert_adj_to_sparse(adj_data, sparse_data):
	for node, neighbors in sorted(adj_data.items()):
		for nid in sorted(neighbors):
			if node < nid:
				sparse_data += [[node, nid, 1]]
			else:
				sparse_data += [[nid, node, 1]]


def write_sparse_file(filename, nodes, edges, sparse_data):
	ff = open(filename, 'w')
	header = "%%MatrixMarket matrix coordinate real general\n%\n"
	head = "%s%d\t%d\t%d\n" % (header, nodes, nodes, edges)
	ff.write(head)
	for node in sparse_data:
		record = "%d\t%d\t%d\n" % (node[0], node[1], node[2])
		ff.write(record)
	ff.close()


infile, outfile = parse_arguments(sys.argv[1:])
adj_data = {}
nnodes, nedges = read_adj_file(infile, adj_data)
print "nodes = %d, edges = %d" % (nnodes, nedges)
sparse_data = []
convert_adj_to_sparse(adj_data, sparse_data)
write_sparse_file(outfile, nnodes, nedges, sparse_data)
