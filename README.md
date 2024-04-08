<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Installing the python library for local development](#installing-the-python-library-for-local-development)
  - [Create a conda environment](#create-a-conda-environment)
  - [Importing functions into your own environment](#importing-functions-into-your-own-environment)
  - [Adding new dependencies to the environment](#adding-new-dependencies-to-the-environment)
- [Tutorial: Using `crocotools_py` to postprocess CROCO model output](#tutorial-using-crocotools_py-to-postprocess-croco-model-output)
- [Tutorial: Run a hindcast simulation locally](#tutorial-run-a-hindcast-simulation-locally)
- [Tutorial: Run a forecast simulation locally](#tutorial-run-a-forecast-simulation-locally)
- [Tutorial: Use docker containers to run a forecast simulation locally](#tutorial-use-docker-containers-to-run-a-forecast-simulation-locally)
- [Tutorial: Set up a server to run the forecast workflow](#tutorial-set-up-a-server-to-run-the-forecast-workflow)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used for our [CROCO](https://www.croco-ocean.org/) related model development.

Directories in the repository: 
- `crocotools_mat`:    some matlab code which is used in conjunction with the official croco-tools to support some of SOMISANA's model development, particularly the pre-processing. 
- `crocotools_py`:     python functions for some postprocessing, validation and plotting CROCO model output
- `download`:          python functions for downloading data used for forcing our CROCO configurations
- `configs`:           configurations used for SOMISANA's hindcast and forecast simulations (see README's in the sub-directories)
- `.github/workflows`: github workflows for running our forecast models operationally on a server set up for this purpose on MIMS (see `run_ops.yml`)

This is largely a redesign of [this](https://github.com/SAEON/somisana) amazing repo (big shout out to Zach Smith), but the new design is more 'model-centric', and we don't include anything website related, which was integrated into the model development in the previous design. 

# Tutorial: installing the python library for local use and/or development

Start off by cloning this repo to your local machine (do this wherever you would like the code):
`git clone git@github.com:SAEON/somisana-croco.git`

Then navigate to the root directory of the repo:
`cd somisana-croco`

## Create a new conda environment

The easiest way to includ the python functions in your own scripts would be to create a new conda environment called 'somisana\_croco' with all the required dependencies already built in. This can be done using the `environment.yml` file in the root directory of the repo (this will ensure that all the required dependencies are added in one step):
```sh
mamba env create -f environment.yml
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

Then activate the environment
```sh
conda activate somisana_croco
```

and install the python library so it's available in your environment:
```sh
pip install --no-deps -e .
```

## Prefer to use your own environment?

If you'd prefer to use your own existing python environment, that's fine too. You just need to make sure you install any dependencies needed by the functions you want to use, and you'll need to do `pip install --no-deps -e .` from the root directory in this repo with your environment activated.


## Adding new dependencies to the environment

If you develop new functions which need dependencies not already in `environment.yml`, please add them to `environment.yml`, and then do:
```sh
mamba env update -f environment.yml --prune
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

The same would apply if others have updated `environment.yml`, and you want to get access to the new functions. You would need to `git pull` to get the latest changes and then run run the update line above. 

# Tutorial: Postprocessing CROCO model output 

## Importing functions into your own environment

Once you've installed the library in your python environment, you should be able to get access to the functions in your own python script like this:
```sh
import crocotools_py.postprocess as post
import crocotools_py.plotting as crocplot
import crocotools_py.validation as val
```

For this tutorial we'll only need to `import crocotools_py.postprocess as post`, and then all functions in `crocotools_py/postprocess.py` are available to use in your own script via `post` e.g. `post.get_var()`. 

## Extracting a variable as an xarray dataarray

The `get_var()` function does most of the heavy lifting, and is called by many other functions. Have a look at the inputs for this function inside `crocotools_py/postprocess.py`. It is basically designed to allow the user to subset the croco output in a few useful ways. The default behaviour is that all the data for a file will be extracted and returned via an xarray dataarray:

```sh
import crocotools_py.postprocess as post
from datetime import datetime

# provide the CROCO files - can be a single file or it can be a file pattern, including wild cards
fname=<your-croco-file(s)>

# define the time origin for the CROCO output files (it's not automatically included in the metadata!)
# important that this corresponds to Yorig used in setting up the CROCO model
ref_date = datetime(1993, 1, 1, 0, 0, 0) # important that this corresponds to Yorig used in setting up the CROCO model
da = post.get_var(fname, 'temp', ref_date = ref_date)
# note the 'time' coordinate contains the real times, which are computed using the ref_date input:
da.time.values
```

extracting all the data can take a long time and can be very large, leading to insufficient memory in the case of big hindcasts so it's useful to rather extract a subset that you may want.

## Extracting a subset in time

This is done using the `tstep` optional input

```sh
# you can spcify a single time-step (zero base index):
da = post.get_var(fname, 'temp', tstep = 45, ref_date = ref_date)
# a range of time-steps:
da = post.get_var(fname, 'temp', tstep = [0,45], ref_date = ref_date)
# a single datetime (the nearest available time-step is found):
da = post.get_var(fname, 'temp', tstep = datetime(1993,2,5), ref_date = ref_date)
# a range of datetimes:
da = post.get_var(fname, 'temp', tstep = [datetime(1993,1,1,12,0,0), datetime(1993,2,15,12,0,0)], ref_date = ref_date)
```

## Extracting a horizontal slice

This is done using the `level` optional input

```sh
# if 'level' >= 0, then a sigma level is extracted. 
# So to extract the bottom layer:
da = post.get_var(fname, 'temp', level = 0, ref_date = ref_date)
# For the surface layer you would need to specify the length of the s_rho dimension minus 1 (due to zero based indexing)
level_surf = np.shape(post.get_var(fname, 'temp', tstep = 0))[0]-1 # many ways to skin this cat
da = post.get_var(fname, 'temp', level = level_surf, ref_date = ref_date)
# if 'level'<0 then it is assumed to be the depth below the water surface (negative downward)
# This requires vertical interpolation to the specified depth, and can take some time to compute for large datasets
# The get data at 100m depth
da = post.get_var(fname, 'temp', level = -100, ref_date = ref_date)
```

## Extracting velocity vector data

The above commands should only be used for extracting data on the 'rho' grid. The u,v components of the horizontal velocity vector are on on their own grids, and are also grid-aligned, not east-north-aligned, so they must be rotated according to the local grid angle. To handle this, we use a dedicated function `get_uv()`

```sh
da_u, da_v = post.get_uv(fname, ref_date = ref_date) # note no 'var' input by definition
# you can apply all the same subsetting as for the `get_var()` function, e.g.
da_u, da_v = post.get_uv(fname, tstep = 45, level = -100, ref_date = ref_date)
```

## Extracting time-series data

The simplest way of doing this is to just specify the `xi_rho`, `eta_rho` indices of the grid cell you want

```sh
da = post.get_var(fname, 'temp', level = level_surf, eta_rho = 50, xi_rho = 50, ref_date = ref_date)
```

But this isn't so useful, particularly for curvilinear grids. We often need to extract data for a specific coordinate e.g.

```sh
lon = 18.4
lat = -34.5
```

To do this we have a dedicated function, called `get_ts()`. E.g to extract the temperature at the provided coordinate and a depth of 100 m:

```sh
ds = post.get_ts(fname, 'temp', lon, lat, ref_date, 
                 time_lims = [datetime(1993,1,1,12,0,0), datetime(1993,2,15,12,0,0)],
                 depths = -100)
```

The output is a dataset which includes both the data and the model depth at the nearest model output location to the specified coordinates. The `depths` input works the same as `level` above and `time_lims` works the same as `tstep` above (something we should probably reconcile!). If the specified depth is lower than the model bottom depth, `nans` are returned, unless you specify `default_to_bottom=True`, then the bottom layer is extracted e.g.

```sh
ds = post.get_ts(fname, 'temp', lon, lat, ref_date, depths = -10000, 
                 default_to_bottom=True)
```

We have a dedicated function for extracting a time-series of the horizontal velocity vectors

```sh
ds = post.get_ts_uv(fname, lon, lat, ref_date, 
                 time_lims = [datetime(1993,1,1,12,0,0), datetime(1993,2,15,12,0,0)],
                 depths = -100)
```

## Extracting profile data

The `get_ts()` and `get_ts_uv()` functions will return data through the water column, depending on the 'depths' input.  If `depths` is not defined, then all sigma levels are extracted and the output dataset returns a `depth` variable which includes the depth of the sigma levels at each time-step:

```sh
ds = post.get_ts(fname, 'temp', lon, lat, ref_date)
```

Alternatively, you can specify your own vertical levels, in which case vertical interpolation is performed

```sh 
depths = [0,-10,-100,-99999]
# if zero is contained in a list of negative numbers, then it is treated as the surface layer
# if -99999 is contained in a list of negative numbers, it is assumed to represent the bottom layer in the model
ds = post.get_ts(fname, 'temp', lon, lat, ref_date, depths = depths)
# Similarly for the velocity vectors
ds = post.get_ts_uv(fname, lon, lat, ref_date, depths = depths)
```

There is also a wrapper function which extracts multiple variables of interest into a single dataset. By default `zeta`, `temp`, `salt`, `u` and `v` are extracted, but you can add variables to the `vars` input if you like:

```sh    
ds = post.get_ts_multivar(fname,lon,lat,ref_date,depths=depths)
```

## Extracting a section

Still to be implemented (any volunteers?)



# Tutorial: Quick plots and animations 

TODO

# Tutorial: validating CROCO output against observations

TODO - Nkululeko

# Tutorial: Running selected python functions via the Command Line Interface (cli.py)

TODO

# Tutorial: Run a hindcast simulation locally

TODO

# Tutorial: Run a forecast simulation locally

TODO

# Tutorial: Introduction to our docker images

TODO

# Tutorial: Use docker images to run a forecast simulation locally

TODO

# Tutorial: Set up a server to run the forecast workflow using Github Actions

TODO
