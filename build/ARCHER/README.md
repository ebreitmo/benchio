# benchio Build: ARCHER

Tested under all three programming environments on ARCHER: Cray, GNU, and Intel.

For an ARCHER build load the required modules:

```
module load cray-netcdf-hdf5parallel/4.4.1.1
module load cray-hdf5-parallel/1.10.0.1
```

Copy ``Makefile.archer`` into ``shared-file/source`` directory

Build:

```bash
cd shared-file/source
make -f Makefile.archer
```

Binary is called ``benchio.x``.
