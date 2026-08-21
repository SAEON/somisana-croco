import crocotools_py.postprocess as post
import xarray as xr
import numpy as np

# we're just going to initialise from rest so the initial file we can create from the grid file and set all params to zero

fname='../../CROCO_GRID/croco_grd.nc'
fname_out='wav.nc'

ds = xr.open_dataset(fname)

zeros=ds.h.values*0

ds = xr.Dataset(
    {
        "WW3_ACHA": (("y", "x"), zeros),
        "WW3__DIR": (("y", "x"), zeros),
        "WW3__OHS":  (("y", "x"), zeros),
        "WW3_T0M1":  (("y", "x"), zeros),
        "WW3_TWOX":  (("y", "x"), zeros),
        "WW3_TWOY":  (("y", "x"), zeros),
        "WW3_TAWX":  (("y", "x"), zeros),
        "WW3_TAWY":  (("y", "x"), zeros),
        "WW3__FOC":  (("y", "x"), zeros),
        "WW3_USSX":  (("y", "x"), zeros),
        "WW3_USSY":  (("y", "x"), zeros),
        "WW3___LM":  (("y", "x"), zeros),
        "WW3__BHD":  (("y", "x"), zeros),
        "WW3_UBRX":  (("y", "x"), zeros),
        "WW3_UBRY":  (("y", "x"), zeros),
    }
)

ds.to_netcdf(fname_out)
