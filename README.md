# benchio
EPCC I/O benchmarking application. Tests write bandwidth to a single shared
file for a given problem size per processor (weak scaling).

Each test is performed 10 times as a default (the number of iterations ```numrep``` can be changed in the input file) and the minimum, maximum and average bandwidth
returned. It is recommended that the maximum bandwidth be considered in most
cases, due to variations in I/O performance from user contention.

Data layout is 3D strided - intended to more closely resemble that of a real
world application than 2D sequential. By default, array size n1xn2xn3 is 256x256x256 (the grid size can be changed in the input file).

Supports POSIX (serial), MPI-IO, HDF5 and NetCDF backends. A run will test all
backends included at compile time.

## Building with selected backends

POSIX (serial), MPI-IO, HDF5 and NetCDF backends can be disabled by commenting
out the relevant `FFLAGS` lines in the `Makefile`. e.g. Commenting out:
`FFLAGS+= -DWITH_NETCDF` will build the application with only MPI-IO and HDF5
support.

## Compilation details

* [ARCHER (EPCC)](build/ARCHER/)
* Athena (HPC Mid+) 
* [Cirrus (EPCC)](build/Cirrus/)
* [CSD3-Skylake (Cambridge)](build/CSD3Skylake/)
* [Thomas (MMM Hub)](build/Thomas/)

For other platforms, edit the `Makefile` to point variables `FC` and `CC` to
MPI compilers (e.g. `mpif90` and `mpicc`), and edit `LFLAGS` to the location of
HDF5 and NetCDF libraries. Then `make clean && make`.

## Job scripts

* [ARCHER (EPCC)](run/ARCHER/)
* Athena (HPC Mid+) 
* [Cirrus (EPCC)](run/Cirrus)
* CSD3-Skylake (Cambridge)
* Thomas (MMM Hub)




## Adjusting volume of data

The array size (per MPI rank) is controlled through the following parameters in
`benchio.F90`:

`integer :: n1 `

`integer :: n2 `

`integer :: n3 `

and are set in `BenchIO-Input.txt`.
