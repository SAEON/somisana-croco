Here is where we prepare all the files to be used in our CROCO configuration and where we run the model

This directory includes a `jobcomp` file for compiling the croco code and a `run_croco_inter.bash` file for initialising the simulation. Both of these scripts use environment variables which are set up in `myenv_inter.bash`.

You'll need to have some directories in place with all of the required inputs, as described below. 

None of the netcdf files are copied to the remote repo, only the scripts used to generate them. 

<pre>
Grid input
----------
`GRID` - this is not intended to be a configurable directory. i.e. if you want a new grid, create a new domain e.g. swcape\_03. The CROCO grid file was generated by the `generate_input.m` script in the `GRID` directory

Compile options
---------------
`C**`
 01 - baseline compile options

Runtime input (\*.in files)
---------------------------
`I**`
 01 - baseline runtime options

Surface and boundary forcing
----------------------------
you'll need to have directory names corresponding to what you have specified in `myenv_inter.bash` i.e. for `ATMOS_BULK` and `OGCM`. For example, in this directory we have a `ERA5` dir and a `GLORYS` dir for surface and boundary forcing, respectively. You can see the matlab scripts called `generate_input.m` in those directories for how the input files can be generated.

</pre>