
# General information about miniconda:

## Overview

Miniconda installers contain the conda package manager and Python.
Once miniconda is installed, you can use the conda command to install any other packages and create environments.

After downloading `Miniconda3-latest-Linux-x86_64.sh` or `Miniconda3-latest-MacOSX-x86_64.sh` you need to run it with: `sh Miniconda3-latest-MacOSX-x86_64.sh`

Add in your .cshrc
```
# Conda environment
source /pathto/.miniconda/etc/profile.d/conda.csh
```
## Main commands:
What version, update conda
```
conda --version
conda update conda
```
Create new environment myenv
```
conda create --name myenv python
```
Switch to another environment
```
conda activate myenv
```
To change your path from the current environment back to the root (or source_deactivate in csh)
```
conda deactivate
```
List all environments
```
conda info --envs
```
Delete an environment
```
conda remove --name myenv --all
```
View a list of packages and versions installed in an environmentSearch for a package
```
conda list
```
Check to see if a package is available for conda to install
```
conda search packagename
```
Install a new package
```
conda install packagename
```
Remove conda
```
rm -rf /pathto/miniconda3 
```
where machine is the name of your computer and username is your username.


## Install a package from Anaconda.org

For packages that are not available using conda install, we can next look on Anaconda.org. Anaconda.org is a package management service for both public and private package repositories. Anaconda.org is a Continuum Analytics product, just like Anaconda and Miniconda.

In a browser, go to http://anaconda.org. We are looking for a package named “pestc4py”
There are more than a dozen copies of petsc4py available on Anaconda.org, first select your platform, then you can sort by number of downloads by clicking the “Downloads” heading.

Select the version that has the most downloads by clicking the package name. This brings you to the Anaconda.org detail page that shows the exact command to use to download it:

Check to see that the package downloaded
```
conda list
```

## Install a package with pip

For packages that are not available from conda or Anaconda.org, we can often install the package with pip (short for “pip installs packages”).

## Exporting environment

```
conda env export > environment.yml on a machine
conda env create -f environment.yml -n $ENV_NAME on the new machine
```


