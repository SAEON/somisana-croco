# scripts to download data for forcing SOMISANA models

In our workflow we would typically download data to a directory called `DATASETS_CROCOTOOLS` in the root of this repo. Doing this will mean that the `crocotools_param.m` files in this repo will be pointing to the right place.

We are just extending the standard external datasets from [croco website](https://www.croco-ocean.org/download-2/) as [DATASETS_CROCOTOOLS.tar.gz](https://data-croco.ifremer.fr/DATASETS/DATASETS_CROCOTOOLS.tar.gz). To do some local development you can start by downwloading this directory and moving it into the root of this repo, and then add new subdirectories such as `GLORYS` or `ERA5` in there. `DATASETS_CROCOTOOLS` is added to .gitignore so nothing you do in there will be tracked in this repo.

To download datasets locally you can set up your own python script to call the download function you want, or you can use the CLI straight in a bash script

