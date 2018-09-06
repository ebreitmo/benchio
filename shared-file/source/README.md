# The Input File

The input file is called BenchIO-Input.txt.
The input is as follows:

````
   # Problem size, dim1:         n1
   # Problem size, dim2:         n2
   # Problem size: dim3:         n3
   # Number of iterations:       numrep
   # Name of output directory:   filedir
   # Overwrite preprocessor values (0 or 1)
                                 WITHSERIAL                 
                                 WITHMPIIO
                                 WITHHDF5
                                 WITHNETCDF
````

The setting of the preprocessor values in the ```Makefile``` for the compiler flags will be overwritten by the entries in the input file. 
A warning will be issued that there is a discrepancy between what is set in the input file and the Makefile.
