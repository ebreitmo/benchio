# The Input File

The input file is called BenchIO-Input.txt.
The input is as follows:

````
   # Problem size, dim1:         n1
   # Problem size, dim2:         n2
   # Problem size: dim3:         n3
   # Number of iterations:       numrep
   # Name of output directory:   filedir
   # Preprocessor values (0 or 1)
                                 WITH_SERIAL                 
                                 WITH_MPIIO
                                 WITH_HDF5
                                 WITH_NETCDF
````

The setting of the preprocessor values will not override the settings in the ```Makefile``` for the compiler flags. 
Only a warning will be issued that there is a discrepancy between what is set in the input file and the Makefile.
