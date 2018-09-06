# Running

For example, to test performance with maximum striping on Lustre:

`mkdir -p benchio_files`

`lfs setstripe -c -1 benchio_files`

`aprun -n <NUMBER_OF_PROCESSORS> ./benchio.x`

Explanation of commands follows:

The application expects a working directory with subdirectory `benchio_files`.
An I/O error will be thrown if this subdirectory is not present.

Under Lustre, the appropriate striping patterns should be set on this
subdirectory using the `lfs setstripe` command.

The `benchio.x` executable should be launched with the platform job launcher,
e.g. `aprun` or `mpirun`. Results are written to STDOUT.

Complete sample run scripts are included in `source/run_scripts`.