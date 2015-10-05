/**
 *  Project:
 *
 *  File: extract_max_level_cell.cpp
 *  Created: Jul 24, 2014
 *  Modified: Thu 24 Jul 2014 09:00:39 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <fstream>
#include <netcdfcpp.h>

#include "netcdf_utils.h"

bool extract_max_level_cell(std::string filename, long int& ncells, int* &max_levels) {
	// open the netcdf grid file
	#ifdef _64BITOFFSET
		NcFile ncid(filename.c_str(), NcFile::ReadOnly, NULL, 0, NcFile::Offset64Bits);
	#else
		NcFile ncid(filename.c_str(), NcFile::ReadOnly);
	#endif

	// read ncells
	NcDim *dim_id = ncid.get_dim("nCells");
	ncells = dim_id->size();
	NcToken cell_dim_token = dim_id->name();

	max_levels = new (std::nothrow) int[ncells];
	if(max_levels == NULL) return false;

	NcVar *var_id = ncid.get_var("maxLevelCell");
	var_id->get(max_levels, ncells);

	ncid.close();

	return true;
} // extract_max_level_cell()


bool write_data(std::string filename, int* data, long int size) {

	std::ofstream outfile(filename);

	for(long int i = 0; i < size; ++ i) outfile << data[i] << std::endl;

	outfile.close();

	return true;
} // write_data()


int main(int narg, char** args) {

	if(narg != 3) {
		std::cout << "usage: extract <ncgridfile> <outfile>" << std::endl;
		return 0;
	} // if

	std::string infile(args[1]);
	std::string outfile(args[2]);

	long int ncells = 0;
	int* max_levels = NULL;
	if(!extract_max_level_cell(infile, ncells, max_levels)) std::cerr << "error!" << std::endl;
	std::cout << "Number of cells = " << ncells << std::endl;
	write_data(outfile, max_levels, ncells);

	delete[] max_levels;

	return 0;
} // main()
