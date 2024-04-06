# scripts to download data for forcing SOMISANA models

For developing hindcast simulations I've been downloading data to a directory called `DATASETS_CROCOTOOLS` in the root of this repo. Doing this will mean that the `crocotools_param.m` files in this repo will be pointing to the right place. `DATASETS_CROCOTOOLS` is added to the `.gitignore` file so nothing you do in there will be tracked in this repo. But of course you are free to download data to wherever you want.

I've just been extending the standard [DATASETS_CROCOTOOLS](https://data-croco.ifremer.fr/DATASETS/DATASETS_CROCOTOOLS.tar.gz) directory which you can download from the [croco website](https://www.croco-ocean.org/download-2/). So to do some local model development you could start by downwloading this directory and moving it into the root of this repo, and then add new subdirectories such as `GLORYS` or `ERA5` in there, and set up your own python scripts to call the download functions in this repo. Again, `DATASETS_CROCOTOOLS` is added to the `.gitignore` file so nothing you do in there will be tracked in this repo.


