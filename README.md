<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [SOMISANA croco related tools](#somisana-croco-related-tools)
  - [Installing the python library for local development](#installing-the-python-library-for-local-development)
    - [Create a conda environment](#create-a-conda-environment)
    - [Importing functions into your own environment](#importing-functions-into-your-own-environment)
  - [Adding new dependencies to the environment](#adding-new-dependencies-to-the-environment)
  - [Developing new configurations](#developing-new-configurations)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# SOMISANA croco related tools

This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used for our [CROCO](https://www.croco-ocean.org/) related model development.

To clone this repo to your local machine:
`git clone git@github.com:SAEON/somisana-croco.git`

Directories in the repository: 
- `crocotools_mat`:    some matlab code which is used in conjunction with the official croco-tools to support some SOMISANA related activities. 
- `crocotools_py`:     python functions for some postprocessing, validation and plotting CROCO model output
- `download`:          python functions for downloading global data used for forcing our CROCO configurations
- `configs`:           configurations used for SOMISANA's hindcast and forecast simulations
- `.github/workflows`: github workflows for running our forecast models operationally (see `run_ops.yml`)

# Installing the python library for local development

## Create a conda environment

The easiest way to include some of the python functions in your own scripts would be to create a new conda environment called 'somisana\_croco' with all the required dependencies already built in. This can be done using the `environment.yml` file:
```sh
mamba env create -f environment.yml
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

Then activate the environment
```sh
conda activate somisana_croco
```

and install the local libraries in this repo into your environment:
```sh
pip install --no-deps -e .
```
Now you do not need to add this directory to your PYTHONPATH as long as you have the somisana\_croco environment activated.

If you'd prefer to use your own python environment, you'll need to install the dependencies needed by the functions you want to use, you could just do `pip install --no-deps -e .` from the root directory in this repo.

## Importing functions into your own environment

If you've got the 'somisana\_croco' environment activated, you should be able to import the libraries you want like this:
```sh
import crocotools_py.postprocess as post
import crocotools_py.plotting as crocplot
import crocotools_py.validation as val
```

## Adding new dependencies to the environment

If you develop new functions which need dependencies not already in `environment.yml`, please add them to `environment.yml`, and then do:
```sh
mamba env update -f environment.yml --prune
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

The same would apply if others have updated `environment.yml`, and you want to get access to the new functions. You would need to `git pull` to get the latest changes and then run run the update line above. 

# Tutorial: Using `crocotools_py` to postprocess CROCO model output 

TODO

# Tutorial: Run a hindcast simulation locally

TODO

# Tutorial: Run a forecast simulation locally

TODO

# Tutorial: Run a dockerised forecast simulation

TODO


