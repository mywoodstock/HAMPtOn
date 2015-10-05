##
 #  Project:
 #
 #  File: convert_patoh_to_mpas.py
 #  Created: Mar 08, 2014
 #  Modified: Thu 03 Apr 2014 06:48:42 AM PDT
 #
 #  Author: Abhinav Sarje <asarje@lbl.gov>
 ##

import sys, getopt

def parse_arguments(argv):
	try:
		opts, args = getopt.getopt(argv, "i:o:")
		if len(opts) < 2:
			raise getopt.GetoptError("Give arguments")
	except getopt.GetoptError:
		print "convert_patoh_to_mpas.py -i <inputfile> -o <outputfile>"
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
	return inputfile, outputfile


def read_patoh_file(filename, data):
	ff = open(filename)
	while True:
		line = ff.readline()
		if not line: break
		parts = map(int, line.strip().split())
		data += parts
	ff.close()


def write_mpas_file(filename, data):
	ff = open(filename, 'w')
	for part in data:
		record = "%d\n" % part
		ff.write(record)
	ff.close()


infile, outfile = parse_arguments(sys.argv[1:])
data = []
read_patoh_file(infile, data)
write_mpas_file(outfile, data)
