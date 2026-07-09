import crocotools_py.postprocess as post
import xarray as xr

fname='../../CROCO_GLORYS/croco_ini_GLORYS_Y2010M01.nc'
time_idx=0 # there is only one time index in the ini file anyway
fname_out='oce.nc'

ds = xr.open_dataset(fname)

ds = ds.isel(time=time_idx)
ds = ds.isel(s_rho=-1)

u=ds.u.values
u=post.u2rho(u)
u=u[1:-1,1:-1]

v=ds.v.values
v=post.v2rho(v)
v=v[1:-1,1:-1]

temp=ds.temp.values
temp=temp[1:-1,1:-1]

ssh=ds.zeta.values
ssh=ssh[1:-1,1:-1]

ds = xr.Dataset(
    {
        "CROCO_EOCE": (("eta_rho", "xi_rho"), u),
        "CROCO_NOCE": (("eta_rho", "xi_rho"), v),
        "CROCO_SST":  (("eta_rho", "xi_rho"), temp),
        "CROCO_SSH":  (("eta_rho", "xi_rho"), ssh),
    }
)

ds.to_netcdf(fname_out)
