# -*- coding: UTF-8 -*-
#

import numpy as np
import xarray as xr
from io_xarray import return_xarray_dataarray, return_xarray_dataset

second2day = 1. / 86400.

# path = "./"
path = "/home/datawork-lops-osi/cmenesg/moz/moz_256x120_ts5_tb0_hc10_new_skpp_bkpp/t1/CROCO_FILES/"
# path = "/home/datawork-lops-osi/cmenesg/moz/moz_512x240_ts5_tb0_hc10_new/t1/CROCO_FILES/"
keymap_files = {
    'coordinate_file': path + "moz_his.nc",
    'metric_file': path + "moz_his.nc",
    'mask_file': path + "moz_his.nc",
    'variable_file': path + "moz_his.nc"
}

keymap_dimensions = {
    'x_rho': 'x_r',
    'y_rho': 'y_r',
    'x_u': 'x_u',
    'y_u': 'y_r',
    'x_v': 'x_r',
    'y_v': 'y_v',
    'x_w': 'x_r',
    'y_w': 'y_r',
    's_rho': 'z_r',
    's_w': 'z_w',
    'time_counter': 't'
}


keymap_coordinates = {
    'nav_lon_rho': 'lon_r',
    'nav_lat_rho': 'lat_r',
    'nav_lon_u': 'lon_u',
    'nav_lat_u': 'lat_u',
    'nav_lon_v': 'lon_v',
    'nav_lat_v': 'lat_v',
    'time_instant': 'time'
}


keymap_variables = {
    'ssh': 'ssh',
    'u': 'u',
    'v': 'v',
    'w': 'w',
    'temp': 'temp',
    'salt': 'salt',
    'rho': 'rho'
}

keymap_metrics = {
    'pm': 'dx_r',
    'pn': 'dy_r',
    'theta_s': 'theta_s',
    'theta_b': 'theta_b',
    'Vtransform': 'scoord',
    'hc': 'hc',
    'h': 'h',
    'f': 'f'
}

keymap_masks = {
    'mask_rho': 'mask_r'
}

# Variables holders for croco


class CrocoWrapper(object):
    """This class create the dictionnary of variables used for creating a
    generic grid from xios croco output files.
    """
    def __init__(self, chunks=None, mask_level=0):

        self.keymap_files = keymap_files
        self.keymap_dimensions = keymap_dimensions
        self.keymap_coordinates = keymap_coordinates
        self.keymap_variables = keymap_variables
        self.keymap_metrics = keymap_metrics
        self.keymap_masks = keymap_masks

        self.chunks = chunks
        self.mask_level = mask_level
        self.coords = {}
        self.metrics = {}
        self.masks = {}

        self.define_coordinates()
        self.define_dimensions(self.dscoord)
        self.define_metrics()
        self.define_masks()
        self.define_variables()
        self.parameters = {}
        self.parameters['chunks'] = chunks

    def _get(self, *args, **kwargs):
        return return_xarray_dataarray(*args, **kwargs)

    def _get_date(self, tindex):
        return self.coords['time'].values[tindex]

    def change_dimensions(self, ds):
        for key, val in self.keymap_dimensions.items():
            try:
                ds = ds.rename({key: val})
            except Exception:
                pass
        return ds

    def change_coords(self, ds):
        for key, val in self.keymap_coordinates.items():
            try:
                ds = ds.rename({key: val})
            except Exception:
                pass
        return ds

    def change_variables(self, ds):
        for key, val in self.keymap_variables.items():
            try:
                ds = ds.rename({key: val})
            except Exception:
                pass
        return ds

    def change_metrics(self, ds):
        for key, val in self.keymap_metrics.items():
            try:
                ds = ds.rename({key: val})
            except Exception:
                pass
        return ds

    def change_mask(self, ds):
        for key, val in self.keymap_masks.items():
            try:
                ds = ds.rename({key: val})
            except Exception:
                pass
        return ds

    def define_dimensions(self, ds):
        self.L = ds.dims['x_r']
        self.M = ds.dims['y_r']
        self.N = ds.dims['z_r']
        self.ntimes = ds.dims['t']

    def define_coordinates(self):
        ds = return_xarray_dataset(self.keymap_files['coordinate_file'])
        ds = self.change_dimensions(ds)
        ds = self.change_coords(ds)
        self.dscoord = ds
        lon_r = self._get(self.dscoord, 'lon_r', chunks=self.chunks, decode_times=False).values
        lat_r = self._get(self.dscoord, 'lat_r', chunks=self.chunks, decode_times=False).values
        self.coords['lon_r'] = lon_r
        self.coords['lat_r'] = lat_r
        self.coords['lon_u'] = 0.5 * (lon_r[:, :-1] + lon_r[:, 1:])
        self.coords['lat_u'] = 0.5 * (lat_r[:, :-1] + lat_r[:, 1:])
        self.coords['lon_v'] = 0.5 * (lon_r[:-1, :] + lon_r[1:, :])
        self.coords['lat_v'] = 0.5 * (lat_r[:-1, :] + lat_r[1:, :])
        self.coords['lon_w'] = lon_r
        self.coords['lat_w'] = lat_r

        # time = time - time_origin
        self.coords['time'] = self._get(self.dscoord, 'time', chunks=self.chunks, decode_times=False)
        self.coords['time'].values = np.array(self.coords['time'], dtype='datetime64[D]') - \
            np.array(self.coords['time'].time_origin, dtype='datetime64[D]')
        self.coords['time'].values = self.coords['time'].values / np.timedelta64(1, 'D')

    def define_metrics(self):
        ds = return_xarray_dataset(self.keymap_files['metric_file'])
        ds = self.change_dimensions(ds)
        ds = self.change_coords(ds)
        ds = self.change_metrics(ds)
        self.dsmetrics = ds
        for key, val in self.keymap_metrics.items():
            self.metrics[val] = self._get(self.dsmetrics, val, chunks=self.chunks)
        # Add missing metrics
        # self.metrics['theta_s'] = xr.DataArray(data=[5.])
        # self.metrics['theta_b'] = xr.DataArray(data=[0.])
        # self.metrics['scoord'] = xr.DataArray(data=[2])
        # rad = np.pi/180
        # Rearth = 6.3708e6
        # dx_r = np.zeros_like(self.coords['lon_r'])
        # dx_r[:,:-1] = np.diff(self.coords['lon_r'],axis=1) * rad * Rearth
        # dx_r[:,-1] = dx_r[:,0]
        # dx_r = 1./dx_r
        # dy_r = np.zeros_like(self.coords['lat_r'])
        # dy_r[:-1,:] = np.diff(self.coords['lat_r'],axis=0) * rad * Rearth
        # dy_r[-1,:] = dy_r[0,:]
        # dy_r = 1./dy_r
        # self.metrics['dx_r'] = xr.DataArray(data=dx_r)
        # self.metrics['dy_r'] = xr.DataArray(data=dy_r)

    def define_masks(self):
        ds = return_xarray_dataset(self.keymap_files['mask_file'])
        ds = self.change_dimensions(ds)
        ds = self.change_coords(ds)
        ds = self.change_mask(ds)
        self.dsmask = ds
        for key, val in self.keymap_masks.items():
            try:
                self.masks[val] = self._get(self.dsmask, val, chunks=self.chunks)
            except Exception:
                mask_rho = np.ones_like(self.coords['lon_r'])
                self.masks[val] = xr.DataArray(data=mask_rho)

    def define_variables(self):
        ds = return_xarray_dataset(self.keymap_files['variable_file'])
        ds = self.change_dimensions(ds)
        ds = self.change_coords(ds)
        ds = self.change_variables(ds)
        # Add ssh as variable
        # ds = ds.assign(ssh = xr.DataArray(data=np.zeros((self.ntimes,self.M,self.L)),dims=('t','y_r','x_r')))
        self.dsvar = ds

    def chunk(self, chunks=None):
        """
        Chunk all the variables.
        Parameters
        ----------
        chunks : dict-like
            dictionnary of sizes of chunk along xarray dimensions.
        """
        for dataname in self.variables:
            data = self.variables[dataname]
            if isinstance(data, xr.DataArray):
                self.variables[dataname] = data.chunk(chunks)

    def _scoord2z(self, point_type, ssh, alpha, beta, lonindex=None, latindex=None):
        """
        scoord2z finds z at either rho or w points (positive up, zero at rest surface)
        h          = array of depths (e.g., from grd file)
        theta_s    = surface focusing parameter
        theta_b    = bottom focusing parameter
        hc         = critical depth
        N          = number of vertical rho-points
        point_type = 'r' or 'w'
        scoord     = 'new2008' :new scoord 2008, 'new2006' : new scoord 2006,
                      or 'old1994' for Song scoord
        ssh       = sea surface height
        message    = set to False if don't want message
        """
        def CSF(sc, theta_s, theta_b):
            '''
            Allows use of theta_b > 0 (July 2009)
            '''
            one64 = np.float64(1)

            if theta_s > 0.:
                csrf = ((one64 - np.cosh(theta_s * sc)) /
                        (np.cosh(theta_s) - one64))
            else:
                csrf = -sc ** 2
            sc1 = csrf + one64
            if theta_b > 0.:
                Cs = ((np.exp(theta_b * sc1) - one64) /
                      (np.exp(theta_b) - one64) - one64)
            else:
                Cs = csrf
            return Cs
        #
        try:
            self.scoord
        except Exception:
            self.scoord = 2
        N = np.float64(self.N)
        try:
            theta_s = self.metrics['theta_s'].values
        except Exception:
            theta_s = self.metrics['theta_s']
        try:
            theta_b = self.metrics['theta_b'].values
        except Exception:
            theta_b = self.metrics['theta_b']
        try:
            hc = self.metrics['hc'].values
        except Exception:
            hc = self.metrics['hc']

        if lonindex is not None:
            h = self.metrics['h'].values[:, lonindex - 1:lonindex + 2]
        elif latindex is not None:
            h = self.metrics['h'].values[latindex - 1:latindex + 2, :]
        else:
            h = self.metrics['h'].values
        scoord = self.metrics['scoord'].values

        sc_w = (np.arange(N + 1, dtype=np.float64) - N) / N
        sc_r = ((np.arange(1, N + 1, dtype=np.float64)) - N - 0.5) / N

        if 'w' in point_type:
            sc = sc_w
            # add a level
            N += 1.
        else:
            sc = sc_r

        z = np.empty((int(N),) + h.shape, dtype=np.float64)
        if scoord == 2:
            Cs = CSF(sc, theta_s, theta_b)
        else:
            try:
                cff1 = 1. / np.sinh(theta_s)
                cff2 = 0.5 / np.tanh(0.5 * theta_s)
            except Exception:
                cff1 = 0.
                cff2 = 0.
            Cs = (1. - theta_b) * cff1 * np.sinh(theta_s * sc) +\
                theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5)

        if scoord == 2:
            hinv = 1. / (h + hc)
            cff = (hc * sc).squeeze()
            cff1 = (Cs).squeeze()
            for k in np.arange(N, dtype=int):
                z[k] = ssh + (ssh + h) * (cff[k] + cff1[k] * h) * hinv
        elif scoord == 1:
            hinv = 1. / h
            cff = (hc * (sc - Cs)).squeeze()
            cff1 = Cs.squeeze()
            cff2 = (sc + 1).squeeze()
            for k in np.arange(N) + 1:
                z0 = cff[k - 1] + cff1[k - 1] * h
                z[k - 1, :] = z0 + ssh * (1. + z0 * hinv)
        else:
            raise Exception("Unknown scoord, should be 1 or 2")
        return z.squeeze(), np.float32(Cs)

    def scoord2z_r(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        Depths at rho point
        '''
        return self._scoord2z('r', ssh=ssh, alpha=alpha, beta=beta,
                              lonindex=lonindex, latindex=latindex)[0]

    def scoord2z_w(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        Depths at rho point
        '''
        return self._scoord2z('w', ssh=ssh, alpha=alpha, beta=beta,
                              lonindex=lonindex, latindex=latindex)[0]

    def scoord2z_u(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        Depths at u point
        '''
        depth = self._scoord2z('r', ssh=ssh, alpha=alpha, beta=beta)[0]
        if lonindex is not None:
            return np.squeeze(0.5 * (depth[:, :, lonindex] + depth[:, :, lonindex - 1]))
        elif latindex is not None:
            return np.squeeze(0.5 * (depth[:, latindex, 1:] + depth[:, latindex, :-1]))
        else:
            return np.squeeze(0.5 * (depth[:, :, :-1] + depth[:, :, 1:]))

    def scoord2z_v(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        Depths at v point
        '''
        depth = self._scoord2z('r', ssh=ssh, alpha=alpha, beta=beta)[0]
        if lonindex is not None:
            return np.squeeze(0.5 * (depth[:, 1:, lonindex] + depth[:, :-1, lonindex]))
        elif latindex is not None:
            return np.squeeze(0.5 * (depth[:, latindex, :] + depth[:, latindex - 1, :]))
        else:
            return np.squeeze(0.5 * (depth[:, :-1, :] + depth[:, 1:, :]))

    def scoord2dz_r(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        """
        dz at rho points, 3d matrix
        """
        dz = self._scoord2z('w', ssh=ssh, alpha=alpha, beta=beta,
                            lonindex=lonindex, latindex=latindex)[0]
        return dz[1:] - dz[:-1]

    def scoord2dz_w(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        """
        dz at rho points, 3d matrix
        """
        dz = self._scoord2z('r', ssh=ssh, alpha=alpha, beta=beta,
                            lonindex=lonindex, latindex=latindex)[0]
        return dz[1:] - dz[:-1]

    def scoord2dz_u(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        dz at u points, 3d matrix
        '''
        dz = self.scoord2z_u(ssh=ssh, alpha=0., beta=1., lonindex=None, latindex=latindex)
        return dz[1:] - dz[:-1]

    def scoord2dz_v(self, ssh=0., alpha=0., beta=1., lonindex=None, latindex=None):
        '''
        dz at v points
        '''
        dz = self.scoord2z_v(ssh=ssh, alpha=0., beta=1., lonindex=lonindex, latindex=latindex)
        return dz[1:] + dz[:-1]


# Run the program

if __name__ == '__main__':
    croco = CrocoWrapper(coordinate_file="moz_his.nc", mask_file="moz_his.nc")
