##
 #  Project:
 #
 #  File: convert_sparse_to_patoh.py
 #  Created: Mar 12, 2014
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 #
 #	Description: Convert sparse graph representation (row, col, val) into hypergraph representation
 #					for input to PATOH.
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print "convert_sparse_to_hyper.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
		else:
			print "wrong argument"
			sys.exit(2)
	return inputfile, outputfile


def read_sparse_file(filename, data):
	ff = open(filename)
	head = ff.readline()
	header = ""
	while head.startswith("%"):
		header += head
		head = ff.readline()
	headwords = head.strip().split()
	if len(headwords) == 2:
		nodes = int(headwords[0])
		edges = int(headwords[1])
	elif len(headwords) == 3:
		nodes = int(headwords[0])
		nodes2 = int(headwords[1])
		edges = int(headwords[2])
		if nodes != nodes2:
			print "error: input matrix is not square"
			ff.close()
			sys.exit(2)
	else:
		print "wrong header in input file"
		ff.close()
		sys.exit(2)
	count = 0
	while True:
		line = ff.readline()
		if not line: break
		count += 1
		words = line.strip().split()
		data += [(int(words[0]), int(words[1]), float(words[2]))]
	if count != edges:
		print "error: mismatch in number of edges and data in input file"
		ff.close()
		sys.exit(2)
	ff.close()
	return nodes, edges


## assuming symmetric data is already present
def convert_sparse_to_hyper(sparse_data, hyper_data):
	pins = 0
	for row, col, val in sparse_data:
		if col not in hyper_data:
			hyper_data[col] = []
		hyper_data[col] += [row]
		pins += 1
	return pins


def write_hyper_file(filename, nnodes, nedges, hyper_data):
	ff = open(filename, 'w')
	head = "%d\t%d\t%d\t%d\n" % (1, nnodes, nnodes, nedges)
	ff.write(head)
	count = 1
	for net, nodes in sorted(hyper_data.items()):
		if count != net:
			print "error: missing data?"
			ff.close()
			sys.exit(2)
		record = ""
		for node in nodes:
			record += "%d\t" % node
		record += "\n"
		ff.write(record)
		count += 1
	ff.close()


infile, outfile = parse_arguments(sys.argv[1:])
sparse_data = []
nnodes, nedges = read_sparse_file(infile, sparse_data)
hyper_data = {}		## { net/hyperedge = [ nodeid, ... ], ... }
npins = convert_sparse_to_hyper(sparse_data, hyper_data)
print "nodes = %d, edges = %d, pins = %d" % (nnodes, nedges, npins)
write_hyper_file(outfile, nnodes, npins, hyper_data)
