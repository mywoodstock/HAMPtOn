/**
 *  Project:
 *
 *  File: extract_cell_coords.cpp
 *  Created: Jul 24, 2014
 *  Modified: Thu 24 Jul 2014 09:00:39 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <fstream>
#include <netcdfcpp.h>

#include "netcdf_utils.h"

bool extract_cell_coords(std::string filename, long int& ncells, float* &coords) {
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

	coords = new (std::nothrow) float[3 * ncells];
	if(coords == NULL) return false;
  float *xcoords, *ycoords, *zcoords;
  xcoords = coords; ycoords = xcoords + ncells; zcoords = ycoords + ncells;

	NcVar *var_x = ncid.get_var("xCell");
	NcVar *var_y = ncid.get_var("yCell");
	NcVar *var_z = ncid.get_var("zCell");
	var_x->get(xcoords, ncells);
	var_y->get(ycoords, ncells);
	var_z->get(zcoords, ncells);

	ncid.close();

	return true;
} // extract_cell_coords()


bool write_coords(std::string filename, float* coords, long int size) {

  float *xcoords, *ycoords, *zcoords;
  xcoords = coords; ycoords = xcoords + size; zcoords = ycoords + size;

	std::ofstream outfile(filename);
	for(long int i = 0; i < size; ++ i)
    outfile << xcoords[i] << "\t"
            << ycoords[i] << "\t"
            << zcoords[i] << std::endl;
	outfile.close();

	return true;
} // write_coords()


int main(int narg, char** args) {

	if(narg != 3) {
		std::cout << "usage: extract <ncgridfile> <outfile>" << std::endl;
		return 0;
	} // if

	std::string infile(args[1]);
	std::string outfile(args[2]);

	long int ncells = 0;
	float* coords = NULL;
	if(!extract_cell_coords(infile, ncells, coords)) std::cerr << "error!" << std::endl;
	std::cout << "Number of cells = " << ncells << std::endl;
	write_coords(outfile, coords, ncells);

	if(coords != NULL) delete[] coords;

	return 0;
} // main()
