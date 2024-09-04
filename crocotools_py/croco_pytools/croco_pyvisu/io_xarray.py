import xarray as xr

# Creation of xarray objects


def return_xarray_dataset(filename, chunks=None, decode_times=True, **kwargs):
    """Return an xarray dataset corresponding to filename.
    Parameters
    ----------
    filename : str
        path to the netcdf file from which to create a xarray dataset
    chunks : dict-like
        dictionnary of sizes of chunk for creating xarray.Dataset.
    Returns
    -------
    ds : xarray.Dataset
    """
    return xr.open_dataset(filename, chunks=chunks, **kwargs)


def return_xarray_mfdataset(filename, chunks=None, **kwargs):
    """Return an xarray dataset corresponding to filename which may include
    wildcards (e.g. file_*.nc).
    Parameters
    ----------
    filename : str
        path to a netcdf file or several netcdf files from which to create a
        xarray dataset
    chunks : dict-like
        dictionnary of sizes of chunk for creating xarray.Dataset.
    Returns
    ------
    ds : xarray.Dataset
    """
    return xr.open_mfdataset(filename, chunks=chunks, **kwargs)


def return_xarray_dataarray(ds, varname, chunks=None, **extra_kwargs):
    """Return a xarray dataarray corresponding to varname in filename.
    Parameters
    ----------
    filename : str
        path to the netcdf file from which to create a xarray.DataArray
    varname : str
        name of the variable from which to create a xarray.DataArray
    chunks : dict-like
        dictionnary of sizes of chunk for creating a xarray.DataArray.
    **extra_kwargs
        not used
    Returns
    -------
    da : xarray.DataArray
    """
    # ds = return_xarray_dataset(filename,chunks=chunks)
    try:
        dataarray = ds[varname]
    except Exception:
        dataarray = ds.attrs[varname]
    for kwargs in extra_kwargs:
        dataarray.attrs[kwargs] = extra_kwargs[kwargs]
    return dataarray
