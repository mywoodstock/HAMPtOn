
# HAMPtOn: README #
HAMPtOn is a Halo-Aware Mesh Partitioning and Ordering tool, typically to be used with MPAS meshes.

## Dependencies
  1. NetCDF
  2. CombBLAS
  3. PaToH

## Installation
  To install the tool, perform the following steps in order.
  1. Make sure the above dependencies are installed.
  2. Define the following environment variables:

    NETCDF_DIR=/path/to/netcdf/install/directory
    COMBBLAS_DIR=/path/to/combblas/directory
    PATOH_BIN=/path/to/patoh/binary/location
    BASE_DIR=/path/to/root/of/this/tool
    TEMP_DIR=/path/to/temporary/space

  `BASE_DIR` and `TEMP_DIR` are optional. Default values for these are the current directory and `/tmp/HAMPtOn`, respectively.
