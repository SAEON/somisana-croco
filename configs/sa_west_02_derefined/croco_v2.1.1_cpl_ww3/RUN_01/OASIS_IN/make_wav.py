import crocotools_py.postprocess as post
import xarray as xr
import numpy as np

# this is just a first attempt to see if I can get the model to start running with this initial condition file
# ultimately we will need to use the previous month's restart to make a new file or rather initialise from rest if not a hot start

fname='/home/gfearon/lustre/somisana-ww3/configs/sa_west_02_derefined/RUN_01/output/ww3.201001.nc'
time_idx=0 # use the first time step in the output file since we are initialising from the start of this month
fname_out='wav.nc'

ds = xr.open_dataset(fname)

ds = ds.isel(time=time_idx)

hs=ds.hs.values
dir=ds.dir.values
dir=(270-dir)*np.pi/180 # conversion geographic to trigonometric and degrees to radians
t0m1=ds.t0m1.values
utwo=ds.utwo.values
vtwo=ds.vtwo.values
utaw=ds.utaw.values
vtaw=ds.vtaw.values
foc=ds.foc.values
uuss=ds.uuss.values
vuss=ds.vuss.values
lm=ds.lm.values
bhd=ds.bhd.values

ds = xr.Dataset(
    {
        "WW3_ACHA": (("y", "x"), hs*0+0.0185),
        "WW3__DIR": (("y", "x"), dir),
        "WW3__OHS":  (("y", "x"), hs),
        "WW3_T0M1":  (("y", "x"), t0m1),
        "WW3_TWOX":  (("y", "x"), utwo),
        "WW3_TWOY":  (("y", "x"), vtwo),
        "WW3_TAWX":  (("y", "x"), utaw),
        "WW3_TAWY":  (("y", "x"), vtaw),
        "WW3__FOC":  (("y", "x"), foc),
        "WW3_USSX":  (("y", "x"), uuss),
        "WW3_USSY":  (("y", "x"), vuss),
        "WW3___LM":  (("y", "x"), lm),
        "WW3__BHD":  (("y", "x"), bhd),
        "WW3_UBRX":  (("y", "x"), uuss),
        "WW3_UBRY":  (("y", "x"), vuss),
    }
)

ds.to_netcdf(fname_out)
