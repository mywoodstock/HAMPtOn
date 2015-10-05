##
 #  Project:
 #
 #  File: symmetrize_sparse.py
 #  Created: Apr 28, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 #
 #	Description: Symmetrize a matrix given in sparse graph format (row, col, val)
 ##

import sys, getopt
from sets import Set

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


def read_sparse_file(filename, data):
	ff = open(filename)
	head = ff.readline()		# this has num nodes and num edges
	while head.startswith("%"): head = ff.readline()
	nodes = int(head.split()[0])
	edges = int(head.split()[2])
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		count += 1
		words = map(int, line.strip().split())
		data.add(tuple(words))
	if count != edges:
		print "error: mismatching number of edges and data in the file: %d, %d" % (edges, count)
		sys.exit(2)
	ff.close()
	return nodes, edges


def symmetrize_sparse(sparse_data):
	more_data = Set()
	for edge in sparse_data:
		more_data.add(tuple([edge[1], edge[0], edge[2]]))
	sparse_data.update(more_data)


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
sparse_data = Set()
nnodes, nedges = read_sparse_file(infile, sparse_data)
print "original nodes = %d, edges = %d" % (nnodes, nedges)
symmetrize_sparse(sparse_data)
newnedges = len(sparse_data)
if newnedges != 2 * nedges: print "warning: newnedges != 2 * nedges. hopefully you know what you are doing"
print "symmetrized nodes = %d, edges = %d" % (nnodes, newnedges)
write_sparse_file(outfile, nnodes, newnedges, sparse_data)
