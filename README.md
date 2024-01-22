# SOMISANA croco related tools

`matlab_tools` contains some matlab code which can be used in conjunction with the official croco-tools to support SOMISANA related activities 
Otherwise, the python code in this directory is intended to useful some some postprocessing, validation and plotting

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

If you want to add more packages to this environment, please add them to `environment.yml`,
and you can then update your environment like this 
```sh
mamba env update -f environment.yml --prune
```
or use `conda` instead of `mamba` if you haven't moved to mamba yet

