# datasets used in SOMISANA croco related model development

This directory would typically be used for downloading external datasets for forcing our configurations. 

For example, you could download the standard external datasets from [croco website](https://www.croco-ocean.org/download-2/) as [DATASETS_CROCOTOOLS.tar.gz](https://data-croco.ifremer.fr/DATASETS/DATASETS_CROCOTOOLS.tar.gz), and move all those datasets into this directory for local development (hence the name of this directory)

You can work in this directory locally, and nothing will get pushed to the remote repo, apart from a few files specified to be included - see ../.gitignore. We're tracking any files which might be useful for others, like the scripts for downloading GLORYS and ERA5 

If there is a new download script you'd like to include in the repo, just specify that you don't want to ignore it in the .gitignore file. 
