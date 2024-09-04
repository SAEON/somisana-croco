XCROCO
====== 




Installation
=============

Install miniconda:
Download Miniconda3 (i.e. for python3) from the [conda website](https://conda.io/miniconda.html) and run:
```
./Miniconda3-latest-Linux-x86_64.sh
```

Download the repository:
```
git clone https://gitlab.inria.fr/croco-ocean/croco_tools.git
cd croco_tools/xcroco
```

Install an appropriate conda-environment:
```
conda env create -n xcroco -f doc/environment.yml
```
or
```
conda update -y conda
conda create -n xcroco -c conda-forge -y python=3.10 
conda activate xcroco
conda install -y -c conda-forge dask dask-jobqueue dask-labextension \
            xarray cf_xarray zarr netcdf4 jupyterlab ipywidgets cartopy \
            geopandas nodejs intake-xarray xgcm numba jupyterhub \
            kerchunk pyamg xrft \
            ffmpeg memory_profiler
            
jupyter labextension install @jupyter-widgets/jupyterlab-manager \
                             @pyviz/jupyterlab_pyviz \
                             jupyter-leaflet
```
Install the xcroco project in an editable mode in the conda environment
```
cd croco-tools/xcroco; pip install -e .
```
see also [conda doc](doc/conda.md)

In order to add the environnement to kernels available to jupyter, you need to run:
```
python -m ipykernel install --user --name xcroco --display-name "xcroco"
```


Tutorials:
=========
You can find tutorials in the directory:
. tuto_xcroco.ipynb : several tools/diagnostics available in the library
. tuto_movie.ipynb : a way to make a movie