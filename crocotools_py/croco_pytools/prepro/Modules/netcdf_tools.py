import netCDF4 as netcdf
import numpy as np

def read_nc(filename,varname, indices="[:]"):
    '''
    Read a variable in a netcdf using netCDF4.Dataset

    Inputs:
      filename          NetCDF filename
      varname           Name of the field in the netcdf
      indices           (Default "[:]") Range of the field in each dimensions

    Outputs:
      var               Retreive var the asked for
    '''
    try:
        with netcdf.Dataset(filename,'r') as nc:
            var = eval(''.join(("nc.variables[varname]", indices)))
    except Exception:
        raise
   #
    if 'float32' in var.dtype.name:
        return var.astype(np.float64)
    else:
        return var


def read_nc_mf(filename,varname,indices="[:]",time_dim='time'):
    '''
    Read a variable in multiple netcdf using netCDF4.MFDataset

    Inputs:
      filename          NetCDF filename , need to be glob.glob(...)
      varname           Name of the field in the netcdf
      indices           (Default "[:]") Range of the field in each dimensions
      time_dim          (Default 'time') Aggregation dimension

    Outputs:
      var               Retreive var the asked for
    '''

    # read multiple netcdf. filename need to be glob.glob(...)
    try:
        try:
            # Load over unlimited dimension
            with netcdf.MFDataset(filename) as nc:
                    var =  eval(''.join(("nc.variables[varname]", indices)))
        except Exception:
            try:
                # Load over time dimension
                with netcdf.MFDataset(filename, aggdim=time_dim) as nc:
                    var =  eval(''.join(("nc.variables[varname]", indices)))
            except Exception:
                print("Vars can not be loaded along %s. Please specify time_dim in read_nc_mf." % time_dim)
    except Exception:
        raise
    if 'float32' in var.dtype.name:
        return var.astype(np.float64)
    else:
        return var

