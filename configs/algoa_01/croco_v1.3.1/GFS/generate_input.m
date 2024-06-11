start
% note that now make_GFS_ocims does not include the step where GFS grb files are reformated to a netcdf file
% this gets handled by a new function called reformat_GFS.m
% make_GFS_OCIMS just interpolates the reformatted data to the blk input file
make_GFS_ocims
