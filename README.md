<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Overview](#overview)
- [Tutorial: installing the python library for local use and/or development](#tutorial-installing-the-python-library-for-local-use-andor-development)
  - [Create a new conda environment](#create-a-new-conda-environment)
  - [Prefer to use your own environment?](#prefer-to-use-your-own-environment)
  - [Adding new dependencies to the environment](#adding-new-dependencies-to-the-environment)
- [Tutorial: Postprocessing CROCO model output](#tutorial-postprocessing-croco-model-output)
  - [Importing functions into your own environment](#importing-functions-into-your-own-environment)
  - [Extracting a variable as an xarray dataarray](#extracting-a-variable-as-an-xarray-dataarray)
  - [Extracting a subset in space](#extracting-a-subset-in-space)
  - [Extracting a subset in time](#extracting-a-subset-in-time)
  - [Extracting a horizontal slice](#extracting-a-horizontal-slice)
  - [Extracting velocity vector data](#extracting-velocity-vector-data)
  - [Extracting time-series data](#extracting-time-series-data)
  - [Extracting profile data](#extracting-profile-data)
  - [Extracting a section](#extracting-a-section)
- [Tutorial: Quick plots and animations](#tutorial-quick-plots-and-animations)
- [Tutorial: validating CROCO output against observations](#tutorial-validating-croco-output-against-observations)
- [Tutorial: Running selected python functions via the Command Line Interface (cli.py)](#tutorial-running-selected-python-functions-via-the-command-line-interface-clipy)
- [Tutorial: Run a hindcast simulation locally](#tutorial-run-a-hindcast-simulation-locally)
- [Tutorial: Run a forecast simulation locally](#tutorial-run-a-forecast-simulation-locally)
  - [download data from global operational models](#download-data-from-global-operational-models)
  - [make CROCO foring files](#make-croco-foring-files)
  - [prepare the runtime options](#prepare-the-runtime-options)
  - [compile the code](#compile-the-code)
  - [run the model](#run-the-model)
- [Tutorial: Introduction to our docker images](#tutorial-introduction-to-our-docker-images)
- [Tutorial: Use docker images to run a forecast simulation locally](#tutorial-use-docker-images-to-run-a-forecast-simulation-locally)
- [Tutorial: Set up a server to run the forecast workflow using Github Actions](#tutorial-set-up-a-server-to-run-the-forecast-workflow-using-github-actions)
  - [set up your user on the server](#set-up-your-user-on-the-server)
  - [install docker](#install-docker)
  - [create a 'somisana' user](#create-a-somisana-user)
  - [configure ~/.bashrc for non-interactive login](#configure-bashrc-for-non-interactive-login)
  - [install github runners](#install-github-runners)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# Overview

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used for our [CROCO](https://www.croco-ocean.org/) related model development.

Directories in the repository: 
- `crocotools_mat`:    some matlab code which is used in conjunction with the official croco-tools to support some of SOMISANA's model development, particularly the pre-processing. 
- `crocotools_py`:     python functions for some postprocessing, validation and plotting CROCO model output
- `download`:          python functions for downloading data used for forcing our CROCO configurations
- `configs`:           configurations used for SOMISANA's hindcast and forecast simulations (see README's in the sub-directories for further details)
- `.github/workflows`: github workflows for running our forecast models operationally on a server set up for this purpose on MIMS (`run_ops.yml` is the `main` workflow)

This repo is largely a redesign of [this one](https://github.com/SAEON/somisana) (big shout out to Zach Smith and Matt Carr for their work on this). The new repo is more 'model-centric', and we don't include anything website related, which was integrated into the model development in the previous design.

# Tutorial: installing the python library for local use and/or development

Start off by cloning this repo to your local machine (do this wherever you would like the code):

`git clone git@github.com:SAEON/somisana-croco.git`

Then navigate to the root directory of the repo:
`cd somisana-croco`

## Create a new conda environment

The easiest way to include the python functions in your own scripts would be to create a new conda environment called 'somisana\_croco' with all the required dependencies already built in. This can be done using the `environment.yml` file in the root directory of the repo (this will ensure that all the required dependencies are added in one step):
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
or use `conda` instead of `mamba` if you haven't moved to mamba yet. This will just update your environment, rather than installing everything from scratch. This is preferred over a `conda install ...` command as the `environment.yml` file keeps a record of all the packages needed to run all python functions.

So if some-one else has updated the `environment.yml` file, and you want to get access to the new functions, you would need to update your local repo with `git pull` and then run the update line above. 

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

# provide the CROCO files - can be a single file or it can be a file pattern, including wild cards. Alternatively `fname` can also be an xarray dataset/dataarray already extracted as part of your script.
fname=<your-croco-file(s)>

# define the time origin for the CROCO output files (it's not automatically included in the metadata!)
# important that this corresponds to Yorig used in setting up the CROCO model
ref_date = datetime(1993, 1, 1, 0, 0, 0) # important that this corresponds to Yorig used in setting up the CROCO model
da = post.get_var(fname, 'temp', ref_date = ref_date)
# note the 'time' coordinate contains the real times, which are computed using the ref_date input:
da.time.values
```

The default behaviour of `get_var()` is to use your CROCO output file to extract the grid information. If your output file does not contain the rho grid variables, `get_var()` has an optional `grdname=<your-grid-file>` input variable, which you can use to specify your croco grid file. 

Extracting all the data can take a long time and can be very large, leading to insufficient memory in the case of big hindcasts so it's useful to rather extract a subset that you may want.

## Extracting a subset in space

This is done using the 'subdomain = [lon0,lon1,lat0,lat1]` optional input

```sh
subdomain=[18,19,-34.5,-34]
da = post.get_var(fname, 'temp', ref_date = ref_date, subdomain=subdomain)
```

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

You can run a forecast simulation locally in a similar way to running a hindcast simulation. For this tutorial, we'll provide the instructions to run a 10 day simulation comprised of a 5 day hindcast component and a 5 day forecast component. In this example, we'll use the `swcape_02` domain and we'll use the CMEMS global ocean forecast model run by Mercator for the boundaries, and GFS for the atmospheric forcing.

## download data from global operational models

First step is to download the forcing files. You can download them locally wherever you want (you'll just need to keep track of your directories in the relevant `crocotools_param.m` files when it comes to generating the CROCO forcing files later).

We have a command line interface (cli.py) for running selected python functions from the command line. This is particularly useful in an operational context. You may need to update your environment so that the cli uses the latest versions of the packages (especially copernicusmarine, which gets updated often as not always backwards compatible):

```sh
mamba env update -f environment.yml --prune
```

Then you should be able to run the download functions:

```sh    
export my_repo_dir='~/code/somisana-croco' # or wherever your local repo is

# data are downloaded over a time period ranging from `run_date - hdays` to `run_date + fdays`
# so define these variables
export run_date=$(date -u +'%Y-%m-%d %H:00:00') # or whatever time you want, in this format (our operational system can be initialised 6 hourly, i.e. %H would be 00, 06, 12 or 18)
export hdays=5 # assuming a 5 day hindcast component
export fdays=5 # assuming a 5 day forecast component

# Mercator download
# It doesn't matter where you download the data locally. In this example we'll download data to `${my_repo_dir}/DATASETS_CROCOTOOLS/`. `DATASETS_CROCOTOOLS` is in the .gitignore file for the repo, so anything you do there can't be pushed to the remote repo.
export mercator_dir=${my_repo_dir}/DATASETS_CROCOTOOLS/MERCATOR # or wherever you want to download mercator data
mkdir ${mercator_dir}
python $my_repo_dir/cli.py download_mercator \
                --usrname <your-copernicus-username> \
                --passwd <your-copernicus-password> \
                --domain 11,23,-39,-25 \ # just enough to cover your domain
                --run_date  "${run_date}"\
                --hdays $hdays \
                --fdays $fdays \
                --outputDir ${mercator_dir}

# GFS download
export gfs_dir=${my_repo_dir}/DATASETS_CROCOTOOLS/GFS # or wherever you want to download GFS data
mkdir ${gfs_dir}
python $my_repo_dir/cli.py download_gfs_atm \
                --domain 11,23,-39,-25 \ # just enough to cover your domain
                --run_date  "${run_date}"\
                --hdays $hdays \
                --fdays $fdays \
                --outputDir ${gfs_dir}
```

The GFS download step above gives us a bunch of `.grb` files - one for each hourly time-step. Don't worry about files which can't be downloaded after the `f120.grb` - after 120 forecast hours (5 days), the data are available as 3 hourly intervals, not hourly, so we don't expect data for `f121.grb`, `f122.grb`,`f124.grb` etc. 

We have a matlab routine to convert these `.grb` files into a single nc file. 
TODO: convert this to python!

```sh
# copy a start file to the GFS download dir, so we can define our paths
cd ${my_repo_dir}/DATASETS_CROCOTOOLS/GFS
cp $my_repo_dir/crocotools_mat/fcst/start_GFS.m . 
# you have to open `start_GFS.m` and change the paths to where your repo's are in your local environment (see instructions in file)
# create a directory for writing the CROCO-friendly nc file
mkdir for_croco
# then you can run the matlab command
# (this will read the gfs.env file created during the download step) 
matlab -nodisplay -nosplash -nodesktop -r "start_GFS; reformat_GFS(2000); exit;" # assuming you're using a Yorig of 2000
``` 

## make CROCO foring files

Let's create a local directory (could be anywhere you want), where we want to run the model

```sh
# create a directory locally 
echo run_dir=~/tmp/croco_forecast_YYYYMMDD_HH/ # obviously change YYYMMDD_HH to whatever date you're running
mkdir -p ${run_dir}
cd ${run_dir}
```

Now we'll copy over configuration files from the repo which we'll use to generate forcing files, compile the code and run the model. The configuration we are using in this example is this one:

```sh
echo config_dir=$my_repo_dir/configs/swcape_02/croco_v1.3.1/ 
```

Have a look at the README in that directory - it explains what all the different directories represent. Let's start by copying over the grid file, which is stored in the repo: 

```sh
cd ${run_dir}
cp -r ${config_dir}/GRID .
```

Then copy over the MERCATOR directory, as that is what we're using to initialise and force our model boundaries.

```sh
# copy over the configuration files for generating MERCATOR forcing
cp -r ${config_dir}/MERCATOR .
cd MERCATOR

# you need to edit start.m to make sure the paths point to where your repo's are in your local environment
# you also need to edit `crocotools_param.m`, specifically `CROCOTOOLS_dir` and `DATADIR`
# `DATADIR` should point to `${my_repo_dir}/DATASETS_CROCOTOOLS` in this example, or wherever you downloaded your data

# then create a directory for storing temporary nc files used in processing global data into croco forcing files
mkdir tmp_for_croco
# get the `run_date` in the format required by the `make_MERCATOR_ocimsi` function
export run_date_formatted=$(date -u +'%Y-%m-%d %H') # whatever date you want - just make sure it is in the right format and that it's the same as the one you used in the download steps above 

# and lastly run the matlab command to create the forcing files
matlab -nodisplay -nosplash -nodesktop -r "start; make_MERCATOR_ocims('${run_date_formatted}',$hdays,1); exit;"
```

And then the GFS forcing files

```sh
# copy over the configuration files for generating MERCATOR forcing
cd ${run_dir}
cp -r ${config_dir}/GFS .
cd GFS

# as per the MERCATOR configuration, you need to edit the same paths in start.m and `crocotools_param.m`
# then the matlab script should create the blk forcing file
matlab -nodisplay -nosplash -nodesktop -r "start; make_GFS_ocims; exit;"
```

## prepare the runtime options

There is a I99 directory which includes the runtime options we are using in our forecasts

```sh
# copy over the runtime input files
cd ${run_dir}
cp -r ${config_dir}/I99 .
cd I99

# compute the time parameters, using the `myenv_in.sh` file for this configuration
source myenv_in.sh
HDAYS=5
FDAYS=5
NDAYS=$((HDAYS + FDAYS))
NUMTIMES=$((NDAYS * 24 * 3600 / DT))
NUMAVG=$((NH_AVG * 3600 / DT))
NUMHIS=$((NH_HIS * 3600 / DT))
NUMAVGSURF=$((NH_AVGSURF * 3600 / DT))
NUMHISSURF=$((NH_HISSURF * 3600 / DT))
NUMRST=$((NH_RST * 3600 / DT)) # frequency of output in restart file to be written
RST_STEP=1 # this is relevant if we want to initialise from a specific time-step in a restart file, since we're initialising from a ini file in this example, we'll use the first and only time-step

# do a sed replacement on the template .in file, and write the output to the scratch dir
sed -e 's/DTNUM/'$DT'/' -e 's/DTFAST/'$DTFAST'/' -e 's/NUMTIMES/'$NUMTIMES'/' -e 's/NUMHISSURF/'$NUMHISSURF'/' -e 's/NUMAVGSURF/'$NUMAVGSURF'/' -e 's/NUMHIS/'$NUMHIS'/' -e 's/NUMAVG/'$NUMAVG'/' -e 's/RST_STEP/'$RST_STEP'/' -e 's/NUMRST/'$NUMRST'/' < croco_fcst.in > croco.in

```

## compile the code

We'll use the C01 compile options for this example. Of course, you're free to try out your own compile options, and just copy and rename the C01 to keep track of your changes: 

```sh
# copy over the files for compiling CROCO
cd ${run_dir}
cp -r ${config_dir}/C01 .
cp ${config_dir}/myenv_frcst.sh .
cp ${config_dir}/jobcomp_frcst.sh .
```

You'll need to edit `myenv_frcst.sh` according to your configuration. Then compile the code, providing the directory with the compile options directory as a command line input (this isn't put in the `myenv_frcst.sh` file as it is dynamically defined elsewhere in our operational workflow):

```sh
./jobcomp_frcst.sh C01
```

## run the model

Now we've got all the inputs needed to run the model, we can create a directory to do this. In our operational workflow we use a directory name which reflects all of the inputs. For example, in this example it would be `C01_I99_MERCATOR_GFS`:

```sh
cd ${run_dir}
mkdir -p C01_I99_MERCATOR_GFS/{scratch,output,postprocess}
cd C01_I99_MERCATOR_GFS/scratch
# get the compiled code
cp ../../C01/croco .
# get the runtime input file
cp ../../I99/croco.in .
# get the grid file
cp ../../GRID/croco_grd.nc .
# get the forcing files (using symbolic links to avoid wasting disk space)
ln -s ../../GFS/croco_blk_GFS_20240527_12.nc croco_blk.nc
ln -s ../../MERCATOR/croco_ini_MERCATOR_20240527_12.nc croco_ini.nc
ln -s ../../MERCATOR/croco_clm_MERCATOR_20240527_12.nc croco_clm.nc
```

We have a simple bash script to run the model - you need to give it the directory name where the model will run (as per `jobcomp_frcst.sh` the command line input is needed as the run directory dynamically defined elsewhere in our operational workflow):

```sh
cd ${run_dir}
cp ${config_dir}/run_croco_frcst.bash .
./run_croco_frcst.bash C01_I99_MERCATOR_GFS
```

If all went according to plan, the model would have run and the output moved to `C01_I99_MERCATOR_GFS/output`. The `C01_I99_MERCATOR_GFS/postprocess` directory is reserved for any analysis you may want to do


# Tutorial: Introduction to our docker images

TODO

# Tutorial: Use docker images to run a forecast simulation locally

TODO

# Tutorial: Set up a server to run the forecast workflow using Github Actions

This is largely taken from [here](https://github.com/SAEON/deployment-platform), but we only use parts of it. Notably we aren’t implementing the docker swarm, because we only deploy the operational models on the server, not any of the web-related stuff. 

## set up your user on the server

When you are able to `ssh <user>@<domain>`, make your life easier by [allowing password-less ssh access](https://www.digitalocean.com/community/tutorials/how-to-configure-ssh-key-based-authentication-on-a-linux-server)

Give your user passwordless-sudo access if you want:
```sh
sudo su

# Optionally set the visudo editor to vim
export EDITOR=vim

visudo # Ensure this line: <name> ALL=(ALL) NOPASSWD:ALL
```

## install docker
The server will need to have docker installed, following the [official instructions](https://docs.docker.com/engine/install/ubuntu/). Then add a user to allow docker to run without having use `sudo` every time:
```sh
sudo usermod -aG docker <user>
```

## create a 'somisana' user
The 'somisana' user will be used for running the models operationally on the server. We also add a 'runners' group (note explicit group ID) which helps with permissions for running our docker images:
```sh
sudo groupadd -g 1999 runners 
sudo adduser --disabled-password --shell /bin/bash --gecos "" somisana
sudo usermod -aG docker somisana
sudo usermod -aG runners somisana
```

## configure ~/.bashrc for non-interactive login

Open `~/.bashrc` for the `somisana` user and delete lines that stop non-interactive users from sourcing this file:

```sh
sudo su
su somisana
vi ~/.bashrc
```

Make sure the following is **_NOT_** present:

```txt
case $- in
    *i*) ;;
      *) return;;
esac
```

## install github runners
To get the repo's github workflows to run on a particular server, you need to set up a new 'self-hosted runner'. From the online github repo (https://github.com/SAEON/somisana-croco), navigate to Settings->Actions->Runners->New self-hosted runner. And follow the instructions.

As per the instructions, the first thing you have to do is create a directory for the runner to be deployed. Each runner needs a new directory, and you can name it whatever you want. But I've followed the convention of `runner<runner-number>`:
```sh
sudo su
su somisana
cd /home/somisana
mkdir runner1
# keep adding directories runner2, runner3 as needed
```

Inside the new directory you just created, follow the instructions to download the installer, unpack it, run ./config.sh with your unique token (you would need a new token for each runner you set up). You will get a few prompts:
- Runner group: default
- Runner name: I've been using the `mims<server-number>-runner<runner-number>` convention e.g. mims1-runner1 for the first runner on the first server on mims
- Then add additional labels where it makes sense. For example, I've been adding the `mims1` label for the first server we were given access to. i.e. `mims<server-number>` convention. As the operations expand, we may want to run models on multiple servers on mims, so the naming convention allows for this. 
- Use the default work folder

You can use the same label for different runner names. If your workflow allows for two jobs to be computed at the same time (i.e. one doesn't need to other to first be completed), then having more than one runner with the same label would allow the jobs to be executed in parallel. If the label is only added to a single runner, then you can't use that label to run jobs in parallel (this can actually be used as a strategy to run jobs in series). So the number of runners to add depends on how many workflow jobs you need to be computed in parallel.

As per the instructions, the last step is to run ./run.sh which connects the runner to github, but now you have an open terminal, and if you close the terminal it closes the connection. So I ended up running it as a screen session, which runs it as a background process:
```sh
Screen -S runner1
./run.sh

```

Now you should see your newly created runner on the github repo (https://github.com/SAEON/somisana-croco) if you navigate to Settings->Actions->Runners. And you can run jobs on this server from a github workflow in this repo by using the `runs_on: mims1`, where `mims1` is the label we provided earlier



