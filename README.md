
# HAMPtOn: README

HAMPtOn is a Halo-Aware Mesh Partitioning and Ordering tool, typically to be used with MPAS meshes.

## Dependencies

1. [**NetCDF**] [netcdf]
2. [**CombBLAS**] [combblas]
3. [**METIS**] [metis]
4. [**PaToH**] [patoh]

  [netcdf]: http://www.unidata.ucar.edu/downloads/netcdf    "NetCDF"
  [combblas]: http://gauss.cs.ucsb.edu/~aydin/CombBLAS      "Combinatorial BLAS"
  [patoh]: http://bmi.osu.edu/umit/software.html            "PaToH"
  [metis]: http://glaros.dtc.umn.edu/gkhome/views/metis     "METIS"

  *Note:* MPI library is a requirement of CombBLAS.
  
  *Note:* Either one of [PaToH][patoh] or [METIS][metis] is required. Pros and cons of each are discussed below.

## Installation

  To install the tool, follow the steps given below.

  1. Make sure the above dependencies are installed.
  2. Define the following environment variables:

        PATOH_BIN     = /path/to/patoh/binary/location
        METIS_BIN     = /path/to/metis/binary/location
        NETCDF_DIR    = /path/to/netcdf/install/directory
        COMBBLAS_DIR  = /path/to/combblas/directory
        BASE_DIR      = /path/to/root/of/this/tool
        TEMP_DIR      = /path/to/temporary/space

  `BASE_DIR` and `TEMP_DIR` are optional. Default values for these are the current directory and `/tmp`, respectively. For example, have a look at the file `vars.sh`.

  3. Use `make` to build and install the tools:

        $ make

  The default compiler is `mpicxx`. To use a different MPI C++ compiler, use:

        $ make MPICXX=<compiler>

  4. On successful compilation, a `bin` directory is created containing four binaries/scripts. You only need to worry about and use `hampton`.


## Usage

 The generated `hampton` script when executed will perform a series of steps in order to generate a halo-aware partitioning and SFC-based data re-ordered mesh files. It takes sevaral arguments:

    hampton <grid_file> <graph_file> <num_parts> [<output_prefix> [<ordering_type> [<num_iterations>]]]
 
Here, the required arguments are:
 1. `<grid_file>`: The input grid/mesh file. These are generally the `.nc` files (in netCDF format).
 2. `<graph_file>`: The input file describing the mesh graph. These are generally named `graph.info` in the MPAS inputs.
 3. `<num_parts>`: The number of partitions to create.

The following are optional arguments, provided to give you more control over the process:

 4. `<output_prefix>`: A string which will be used as prefix in the generated output files. Default value is same as the `<graph_file>`.
 5. `<ordering_type>`: A string from the set `{ hilbert, morton, peano, random, xyz, none }`. This specifies the type of data re-ordering to perform. Default value is `hilbert`. A value of `none` will disable ordering.
 6. `<num_iterations>`: The number of iterations to perform in the monte-carlo process for convergence. Default value is 20.

For example:

     $ ./bin/hampton grid.nc graph.info 64 newmesh

  This will generate a output mesh file named `newmesh.64.nc`, a graph file named `newmesh.64.info`, and a partitioning file named `newmesh.64.info.part.64`.

  The generated mesh, graph and partitioning files can then be directly used in your application. For MPAS-Ocean, just update the file names in the `namelist.ocean` and `streams.ocean` input files, and you are good to go.


## PaToH or METIS?
 
 Discussion about pros and cons of both partitioners will come here some day.

