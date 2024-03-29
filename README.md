# SOMISANA croco related tools

`crocotools_mat` contains some matlab code which can be used in conjunction with the official croco-tools to support SOMISANA related activities. 
`crocotools_py` contains some additional python code which is hopefully useful for some postprocessing, validation and plotting.
The files can be easily imported into your code as packages e.g.
```sh
import crocotools_py.postprocess as post
import crocotools_py.plotting as crocplot
import crocotools_py.validation as val
```
... but only if you have an environment set up as described below

## Create the environment

It's easiest to to create a new conda environment called 'somisana\_croco' with all the required dependencies already built in. This can be done using the `environment.yml` file:
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

## Adding new packages to the environment

If you want to add more packages to this environment, please add them to `environment.yml` (pls push commits to this repo),
and you can then update your environment like this 

You'll also need to run this if packages need to be updated, or of course if some-one else has added to the `environment.yml` file

```sh
mamba env update -f environment.yml --prune
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

## Developing new configurations

We are also using this repo for the development of new croco configurations, in the `configs` directory.

