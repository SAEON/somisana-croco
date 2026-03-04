import numpy as np
import Cgrid_transformation_tools as grd_tools
import sigmagrid_tools as sig_tools
import netcdf_tools as nc_tools
import netCDF4 as netcdf
from datetime import datetime
from collections import OrderedDict

class CROCO_grd(object):

    def __init__(self, filename, sigma_params=None):
        print('Reading CROCO grid: %s' %filename )
        self.grid_file = filename
        nc=netcdf.Dataset(filename,'r')
        self.lon = nc_tools.read_nc(filename,'lon_rho')
        self.lat = nc_tools.read_nc(filename,'lat_rho')
        self.lonu =nc_tools.read_nc(filename,'lon_u')
        self.latu = nc_tools.read_nc(filename,'lat_u')
        self.lonv = nc_tools.read_nc(filename,'lon_v')
        self.latv = nc_tools.read_nc(filename,'lat_v')
        self.pm  = nc_tools.read_nc(filename,'pm')
        self.pn  = nc_tools.read_nc(filename,'pn')
        self.maskr = nc_tools.read_nc(filename,'mask_rho')
        self.angle = nc_tools.read_nc(filename,'angle')
        self.h = nc_tools.read_nc(filename,'h')
        self.hraw = nc_tools.read_nc(filename,'hraw')
        self.f = nc_tools.read_nc(filename,'f')
        self.umask= self.maskr[:,0:-1]*self.maskr[:,1:]
        self.vmask= self.maskr[0:-1,:]*self.maskr[1:,:]
        self.pmask= self.umask[0:-1,:]*self.umask[1:,:]
        if sigma_params is not None:
            self.theta_s = np.double(sigma_params['theta_s'])
            self.theta_b = np.double(sigma_params['theta_b'])
            self.hc = np.double(sigma_params['hc'])
            self.N = np.double(sigma_params['N'])
        self.sc_r = None
        nc.close
    def mask3d(self):
        return np.tile(self.maskr, (int(self.N), 1, 1))

    def umask3d(self):
        return np.tile(self.umask, (int(self.N), 1, 1))

    def vmask3d(self):
        return np.tile(self.vmask, (int(self.N), 1, 1))   

    def lonmin(self):
        return np.min(self.lon)

    def lonmax(self):
        return np.max(self.lon)

    def latmin(self):
        return np.min(self.lat)

    def latmax(self):
        return np.max(self.lat)

    def scoord2z_r(self, zeta=0., bdy="", scoord='new2008'):
        '''
        Depths at vertical rho points
        '''
        return sig_tools.scoord2z('r', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[0]


    def Cs_r(self, zeta=0., bdy="", scoord='new2008'):
        '''
        S-coordinate stretching curves at rho points
        '''
        return sig_tools.scoord2z('r', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[1]


    def scoord2z_w(self, zeta=0., bdy="", scoord='new2008'):
        '''
        Depths at vertical w points
        '''
        return sig_tools.scoord2z('w', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[0]

    def Cs_w(self, zeta=0., bdy="", scoord='new2008'):
        '''
        S-coordinate stretching curves at w points
        '''
        return sig_tools.scoord2z('w', zeta=zeta, topo=eval(''.join(("self.h",bdy))), theta_s=self.theta_s, theta_b=self.theta_b,\
                N=self.N,hc=self.hc,scoord=scoord)[1]
    def s_rho(self):
        return ((np.arange(1,self.N+1,dtype=np.float64))-self.N-0.5)/self.N
    def s_w(self):
        return (np.arange(self.N+1,dtype=np.float64)-self.N)/self.N

    def WEST_grid(self,indices="[:,0:2]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_west=eval(''.join(('self.h',indices)))
        self.lon_west=eval(''.join(('self.lon',indices)))
        self.lat_west=eval(''.join(('self.lat',indices)))
        self.lonu_west=eval(''.join(('self.lonu',indices)))
        self.latu_west=eval(''.join(('self.latu',indices)))
        self.lonv_west=eval(''.join(('self.lonv',indices)))
        self.latv_west=eval(''.join(('self.latv',indices)))
        self.maskr_west=eval(''.join(('self.maskr',indices)))
        self.umask_west=eval(''.join(('self.umask',indices)))
        self.vmask_west=eval(''.join(('self.vmask',indices)))
        self.angle_west=eval(''.join(('self.angle',indices)))

        return self

    def EAST_grid(self,indices="[:,-2:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_east=eval(''.join(('self.h',indices)))
        self.lon_east=eval(''.join(('self.lon',indices)))
        self.lat_east=eval(''.join(('self.lat',indices)))
        self.lonu_east=eval(''.join(('self.lonu',indices)))
        self.latu_east=eval(''.join(('self.latu',indices)))
        self.lonv_east=eval(''.join(('self.lonv',indices)))
        self.latv_east=eval(''.join(('self.latv',indices)))
        self.maskr_east=eval(''.join(('self.maskr',indices)))
        self.umask_east=eval(''.join(('self.umask',indices)))
        self.vmask_east=eval(''.join(('self.vmask',indices)))
        self.angle_east=eval(''.join(('self.angle',indices)))

        return self

    def SOUTH_grid(self,indices="[0:2,:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_south=eval(''.join(('self.h',indices)))
        self.lon_south=eval(''.join(('self.lon',indices)))
        self.lat_south=eval(''.join(('self.lat',indices)))
        self.lonu_south=eval(''.join(('self.lonu',indices)))
        self.latu_south=eval(''.join(('self.latu',indices)))
        self.lonv_south=eval(''.join(('self.lonv',indices)))
        self.latv_south=eval(''.join(('self.latv',indices)))
        self.maskr_south=eval(''.join(('self.maskr',indices)))
        self.umask_south=eval(''.join(('self.umask',indices)))
        self.vmask_south=eval(''.join(('self.vmask',indices)))
        self.angle_south=eval(''.join(('self.angle',indices)))

        return self

    def NORTH_grid(self,indices="[-2:,:]"):
        '''
        Defines some vars for Western grid
        '''
        self.h_north=eval(''.join(('self.h',indices)))
        self.lon_north=eval(''.join(('self.lon',indices)))
        self.lat_north=eval(''.join(('self.lat',indices)))
        self.lonu_north=eval(''.join(('self.lonu',indices)))
        self.latu_north=eval(''.join(('self.latu',indices)))
        self.lonv_north=eval(''.join(('self.lonv',indices)))
        self.latv_north=eval(''.join(('self.latv',indices)))
        self.maskr_north=eval(''.join(('self.maskr',indices)))
        self.umask_north=eval(''.join(('self.umask',indices)))
        self.vmask_north=eval(''.join(('self.vmask',indices)))
        self.angle_north=eval(''.join(('self.angle',indices)))

        return self

class CROCO():

    def create_ini_nc(self, filename, grdobj, created_by='make_ini.py',tracers=['temp','salt']):#fillval
        # Global attributes
        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type  = 'CROCO initial file produced by %s' % created_by
        nc.grd_file = grdobj.grid_file
        nc.hc = grdobj.hc
        nc.theta_s = grdobj.theta_s
        nc.theta_b = grdobj.theta_b
        nc.Tcline = grdobj.hc
        nc.Cs_r = grdobj.Cs_r()
        nc.Cs_w = grdobj.Cs_w()
        nc.VertCoordType = 'NEW'

        # Dimensions
        nc.createDimension('xi_rho', grdobj.lon.shape[1])
        nc.createDimension('xi_u', grdobj.lon.shape[1] - 1)
        nc.createDimension('eta_rho', grdobj.lon.shape[0])
        nc.createDimension('eta_v', grdobj.lon.shape[0] - 1)
        nc.createDimension('s_rho', grdobj.N)
        nc.createDimension('s_w', grdobj.N + 1)
        nc.createDimension('time', None)
        nc.createDimension('one', 1)

        # Create the variables and write...
        nc.createVariable('theta_s', 'f', ('one'), zlib=True)
        nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
        nc.variables['theta_s'].units = 'nondimensional'
        nc.variables['theta_s'][:] = grdobj.theta_s

        nc.createVariable('theta_b', 'f', ('one'), zlib=True)
        nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
        nc.variables['theta_b'].units = 'nondimensional'
        nc.variables['theta_b'][:] = grdobj.theta_b

        nc.createVariable('Tcline', 'f', ('one'), zlib=True)
        nc.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
        nc.variables['Tcline'].units = 'meters'
        nc.variables['Tcline'][:] = grdobj.hc

        nc.createVariable('hc', 'f', ('one'), zlib=True)
        nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
        nc.variables['hc'].units = 'meters'
        nc.variables['hc'][:] = grdobj.hc

        nc.createVariable('s_rho', 'f8', ('s_rho'))
        nc.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
        nc.variables['s_rho'].units = 'nondimensional'
        nc.variables['s_rho'].valid_min = -1.
        nc.variables['s_rho'].valid_max = 0.
        nc.variables['s_rho'][:] = grdobj.s_rho()

        nc.createVariable('Cs_rho', 'f8', ('s_rho'), zlib=True)
        nc.variables['Cs_rho'].long_name = 'S-coordinate stretching curves at RHO-points'
        nc.variables['Cs_rho'].units = 'nondimensional'
        nc.variables['Cs_rho'].valid_min = -1.
        nc.variables['Cs_rho'].valid_max = 0.
        nc.variables['Cs_rho'][:] = grdobj.Cs_r()

        nc.createVariable('Cs_w', 'f8', ('s_w'), zlib=True)
        nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at w-points'
        nc.variables['Cs_w'].units = 'nondimensional'
        nc.variables['Cs_w'].valid_min = -1.
        nc.variables['Cs_w'].valid_max = 0.
        nc.variables['Cs_w'][:] = grdobj.Cs_w()

        nc.createVariable('ocean_time', 'f8', ('time'), zlib=True)
        nc.variables['ocean_time'].long_name = 'time since initialization'
        nc.variables['ocean_time'].units     = 'seconds'

        nc.createVariable('scrum_time', 'f8', ('time'), zlib=True)
        nc.variables['scrum_time'].long_name = 'time since initialization'
        nc.variables['scrum_time'].units     = 'seconds'

        nc.createVariable('tstart', 'f8', ('one'), zlib=True)
        nc.variables['tstart'].long_name = 'start processing day'
        nc.variables['tstart'].units     = 'days'

        nc.createVariable('tend', 'f8', ('one'), zlib=True)
        nc.variables['tend'].long_name = 'end processing day'
        nc.variables['tend'].units     = 'days'

        # dictionary for the prognostic variables
        prog_vars = OrderedDict()
        prog_vars['u']    = ['u3d',
                             'initial u-momentum component',
                             'meters second-1']
        prog_vars['v']    = ['v3d',
                             'initial v-momentum component',
                             'meters second-1']
        prog_vars['ubar'] = ['u2d',
                             'initial vertically integrated u-momentum component',
                             'meters second-1']
        prog_vars['vbar'] = ['v2d',
                             'initial vertically integrated v-momentum component',
                             'meters second-1']
        prog_vars['zeta'] = ['rho2d',
                             'initial sea surface height',
                             'meters']
        for trc in tracers:
            if trc == 'temp':
                prog_vars['temp'] = ['rho3d',
                                    'initial potential temperature',
                                    'Celsius']
            elif trc == 'salt':
                prog_vars['salt'] = ['rho3d',
                                     'initial salinity',
                                     'psu']
            else:
                prog_vars[trc] = ['rho3d','','']                

        for varname, value in zip(prog_vars.keys(), prog_vars.values()):

            if 'rho3d' in value[0]:
                dims = ('time', 's_rho', 'eta_rho', 'xi_rho')

            elif 'u3d' in value[0]:
                dims = ('time', 's_rho', 'eta_rho', 'xi_u')

            elif 'v3d' in value[0]:
                dims = ('time', 's_rho', 'eta_v', 'xi_rho')

            elif 'u2d' in value[0]:
                dims = ('time', 'eta_rho', 'xi_u')

            elif 'v2d' in value[0]:
                dims = ('time', 'eta_v', 'xi_rho')

            elif 'rho2d' in value[0]:
                dims = ('time', 'eta_rho', 'xi_rho')

            else: error

            nc.createVariable(varname, 'f8', dims, zlib=True)#fill_value=fillval
            nc.variables[varname].long_name = value[1]
            nc.variables[varname].units     = value[2]

        nc.close()

    def create_bry_nc(self,filename, grdobj, obc_dict, cycle, created_by='make_bry.py',time=None,tracers=['temp','salt']):#fillval
        # Global attributes
        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type = 'CROCO boundary file produced by %s' %created_by
        nc.grd_file = grdobj.grid_file
        nc.hc = grdobj.hc
        nc.theta_s = grdobj.theta_s
        nc.theta_b = grdobj.theta_b
        nc.Tcline = grdobj.hc
        nc.Cs_r = grdobj.Cs_r()
        nc.Cs_w = grdobj.Cs_w()
        nc.VertCoordType = 'NEW'

        # Dimensions
        nc.createDimension('xi_rho', grdobj.lon.shape[1])
        nc.createDimension('xi_u', grdobj.lon.shape[1]-1)
        nc.createDimension('eta_rho', grdobj.lon.shape[0])
        nc.createDimension('eta_v', grdobj.lon.shape[0]-1)
        nc.createDimension('s_rho', grdobj.N)
        nc.createDimension('s_w', grdobj.N+1)
        nc.createDimension('bry_time', None)
        nc.createDimension('one', 1)

        # Create the variables and write...
        nc.createVariable('theta_s', 'f', ('one'))
        nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
        nc.variables['theta_s'].units = 'nondimensional'
        nc.variables['theta_s'][:] = grdobj.theta_s

        nc.createVariable('theta_b', 'f', ('one'))
        nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
        nc.variables['theta_b'].units = 'nondimensional'
        nc.variables['theta_b'][:] = grdobj.theta_b

        nc.createVariable('Tcline', 'f', ('one'))
        nc.variables['Tcline'].long_name = 'S-coordinate surface/bottom layer width'
        nc.variables['Tcline'].units = 'meters'
        nc.variables['Tcline'][:] = grdobj.hc

        nc.createVariable('hc', 'f', ('one'))
        nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
        nc.variables['hc'].units = 'meters'
        nc.variables['hc'][:] = grdobj.hc

        nc.createVariable('s_rho', 'f8', ('s_rho'))
        nc.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
        nc.variables['s_rho'].units = 'nondimensional'
        nc.variables['s_rho'].valid_min = -1.
        nc.variables['s_rho'].valid_max = 0.
        nc.variables['s_rho'][:] = grdobj.s_rho()

        nc.createVariable('Cs_rho', 'f8', ('s_rho'))
        nc.variables['Cs_rho'].long_name = 'S-coordinate stretching curves at RHO-points'
        nc.variables['Cs_rho'].units = 'nondimensional'
        nc.variables['Cs_rho'].valid_min = -1.
        nc.variables['Cs_rho'].valid_max = 0.
        nc.variables['Cs_rho'][:] = grdobj.Cs_r()

        nc.createVariable('Cs_w', 'f8', ('s_w'))
        nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at w-points'
        nc.variables['Cs_w'].units = 'nondimensional'
        nc.variables['Cs_w'].valid_min = -1.
        nc.variables['Cs_w'].valid_max = 0.
        nc.variables['Cs_w'][:] = grdobj.Cs_w()

        nc.createVariable('bry_time', 'f8', ('bry_time'), zlib=True)
        nc.variables['bry_time'].long_name = 'time for boundary data'
        nc.variables['bry_time'].units = 'days'

        # dictionary for the prognostic variables
        prog_vars = OrderedDict()
        prog_vars['u_']    = ['u2',
                              ' boundary u-momentum component',
                              'meters second-1']
        prog_vars['v_']    = ['v2',
                              ' boundary v-momentum component',
                              'meters second-1']
        prog_vars['ubar_'] = ['u1',
                              ' boundary vertically integrated u-momentum component',
                              'meters second-1']
        prog_vars['vbar_'] = ['v1',
                              ' boundary vertically integrated v-momentum component',
                              'meters second-1']
        prog_vars['zeta_'] = ['rho1',
                              ' boundary sea surface height',
                              'meters']
        for trc in tracers:
            if trc == 'temp':
                prog_vars['temp_'] = ['rho2',
                                      ' boundary potential temperature',
                                      'Celsius']
            elif trc == 'salt':
                prog_vars['salt_'] = ['rho2',
                                      ' boundary salinity',
                                      'psu']
            else:
                prog_vars[f"{trc}_"] = ['rho2','','']

        # Loop over boundary
        for boundary, flag in zip(obc_dict.keys(), obc_dict.values()):
            if flag:
                varlabel = '%sern'   % boundary
                for key, value in zip(prog_vars.keys(), prog_vars.values()):
                    if 'rho2' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 's_rho', 'eta_rho')
                        elif boundary=='south' or boundary=='north':
                            dims = ('bry_time', 's_rho', 'xi_rho')
                    elif 'u2' in value[0]:
                        if boundary=='south' or boundary=='north':
                            dims = ('bry_time', 's_rho', 'xi_u')
                        else:
                            dims = ('bry_time', 's_rho', 'eta_rho')
                    elif 'v2' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 's_rho', 'eta_v')
                        else:
                            dims = ('bry_time', 's_rho', 'xi_rho')
                    elif 'u1' in value[0]:
                        if boundary=='south' or boundary=='north':
                            dims = ('bry_time', 'xi_u')
                        else:
                            dims = ('bry_time', 'eta_rho')
                    elif 'v1' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 'eta_v')
                        else:
                            dims = ('bry_time', 'xi_rho')
                    elif 'rho1' in value[0]:
                        if boundary=='east' or boundary=='west':
                            dims = ('bry_time', 'eta_rho')
                        elif boundary=='south' or boundary=='north':
                            dims = ('bry_time', 'xi_rho')
                    else: error

                    varname  = ''.join((key, '%s' % boundary))
                    nc.createVariable(varname, 'f8', dims, zlib=True)
                    nc.variables[varname].long_name = varlabel + value[1]
                    nc.variables[varname].units     = value[2]

        nc.close()

    def create_clim_nc(self, filename, grdobj, cycle, created_by='make_clim.py', tracers=['temp', 'salt']):
    
        # Generic functions for making the various variables
        def create_scalar(nc, name, dtype='f8', dims=('one',), **attrs):
            var = nc.createVariable(name, dtype, dims)
            for k, v in attrs.items():
                setattr(var, k, v)
            return var
        
        def create_coord(nc, name, dims, data=None, **attrs):
            var = nc.createVariable(name, 'f8', dims)
            for k, v in attrs.items():
                setattr(var, k, v)
            if data is not None:
                var[:] = data
            return var
    
        def create_time_variable(nc, name, long_name):
            return create_field(
                nc,
                name,
                (name,),
                long_name=long_name,
                units = 'day', 
                calendar = '360_day',
                cycle_length=0.
            )

        def create_field(nc, name, dims, units, long_name,
                         time_name=None, **extra_attrs):
            var = nc.createVariable(
                name, 'f8', dims,
                zlib=True, complevel=4, shuffle=True
            )
            var.long_name = long_name
            var.units = units
            if time_name:
                var.time = time_name
            for k, v in extra_attrs.items():
                setattr(var, k, v)
            return var
    
        # Create the file
        with netcdf.Dataset(filename, 'w', format='NETCDF4') as nc:
    
            # Global attributes
            nc.date = datetime.now().strftime("%d-%b-%Y")
            nc.clim_file = filename
            nc.grd_file = grdobj.grid_file
            nc.type = f'CROCO climatology file produced by {created_by}'
            nc.hc = grdobj.hc
            nc.theta_s = grdobj.theta_s
            nc.theta_b = grdobj.theta_b
            nc.VertCoordType = 'NEW'
    
            # Dimensions
            eta_rho, xi_rho = grdobj.lon.shape
    
            dims = {
                'xi_u': xi_rho - 1,
                'xi_v': xi_rho,
                'xi_rho': xi_rho,
                'eta_u': eta_rho,
                'eta_v': eta_rho - 1,
                'eta_rho': eta_rho,
                's_rho': grdobj.N,
                's_w': grdobj.N + 1,
                'one': 1,
                'tracer': np.size(tracers)
            }
    
            for name, size in dims.items():
                nc.createDimension(name, size)
    
            time_dims = [
                'tclm_time', 'temp_time',
                'sclm_time', 'salt_time',
                'uclm_time', 'vclm_time',
                'v2d_time', 'v3d_time',
                'ssh_time', 'zeta_time'
            ]
    
            for tdim in time_dims:
                nc.createDimension(tdim, None)
    
            # Classic scalar variables 
            create_scalar(
                nc, 'spherical', 'S1',
                long_name='Grid type logical switch',
                option_T='spherical Cartesian'
            )[:] = 'T'
    
            create_scalar(
                nc, 'Vtransform', 'f8',
                long_name='vertical terrain-following transformation equation'
            )
    
            create_scalar(
                nc, 'Vstretching', 'f8',
                long_name='vertical terrain-following stretching function'
            )
    
            create_scalar(
                nc, 'tstart', 'f8',
                long_name='start processing day',
                units='day'
            )
    
            create_scalar(
                nc, 'tend', 'f8',
                long_name='end processing day',
                units='day'
            )
    
            create_scalar(
                nc, 'theta_s', 'f8',
                long_name='S-coordinate surface control parameter',
                units='nondimensional'
            )
    
            create_scalar(
                nc, 'theta_b', 'f8',
                long_name='S-coordinate bottom control parameter',
                units='nondimensional'
            )
    
            create_scalar(
                nc, 'Tcline', 'f8',
                long_name='S-coordinate surface/bottom layer width',
                units='meter'
            )[:] = grdobj.hc
            
            create_scalar(
                nc, 'hc', 'f8',
                long_name='S-coordinate parameter, critical depth',
                units='meter'
            )
            
            # Vertical coordinates
            create_coord(
                nc, 'sc_r', ('s_rho',),
                data = grdobj.s_rho(),
                long_name='S-coordinate at RHO-points',
                valid_min=-1.,
                valid_max=0.,
                positive='up',
                standard_name='ocean_s_coordinate_g2',
                formula_terms='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
            )
    
            create_coord(
                nc, 'sc_w', ('s_w',),
                data = grdobj.s_w(),
                long_name='S-coordinate at W-points',
                valid_min=-1.,
                valid_max=0.,
                positive='up',
                standard_name='ocean_s_coordinate_g2',
                formula_terms='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
            )
            
            create_coord(
                nc, 'Cs_r', ('s_rho',),
                data = grdobj.Cs_r(),
                long_name='S-coordinate stretching curves at RHO-points',
                units='nondimensional',
                valid_min=-1.,
                valid_max=0.
            )
            
            create_coord(
                nc, 'Cs_w', ('s_w',),
                data = grdobj.Cs_w(),
                long_name='S-coordinate stretching curves at W-points',
                units='nondimensional',
                valid_min=-1.,
                valid_max=0.
            )
            
            # Time variables 
            time_metadata = {
                'tclm_time': 'time for temperature climatology',
                'temp_time': 'time for temperature climatology',
                'sclm_time': 'time for salinity climatology',
                'salt_time': 'time for salinity climatology',
                'uclm_time': 'time for u-momentum climatology',
                'vclm_time': 'time for v-momentum climatology',
                'v2d_time': 'time for 2D velocity climatology',
                'v3d_time': 'time for 3D velocity climatology',
                'ssh_time': 'time for sea surface height',
                'zeta_time': 'time for sea surface height',
            }
    
            for name, long_name in time_metadata.items():
                create_time_variable(nc, name, long_name)
    
            # Tracers 
            tracer_metadata = {
                'temp': ('tclm_time', 'potential temperature', 'Celsius'),
                'salt': ('sclm_time', 'salinity', 'PSU'),
            }
    
            for tracer in tracers:
                if tracer not in tracer_metadata:
                    raise ValueError(f"Unsupported tracer: {tracer}")
                time_name, long_name, units = tracer_metadata[tracer]
    
                create_field(
                    nc, tracer,
                    (time_name, 's_rho', 'eta_rho', 'xi_rho'),
                    units,
                    long_name,
                    time_name=time_name,
                    coordinates=f"lon_rho lat_rho s_rho {time_name}"
                )
    
            # Velocity fields 
            create_field(
                nc, 'u',
                ('uclm_time', 's_rho', 'eta_u', 'xi_u'),
                'meter second-1',
                'u-momentum component',
                time_name='uclm_time'
            )
    
            create_field(
                nc, 'v',
                ('vclm_time', 's_rho', 'eta_v', 'xi_v'),
                'meter second-1',
                'v-momentum component',
                time_name='vclm_time'
            )
    
            create_field(
                nc, 'ubar',
                ('uclm_time', 'eta_u', 'xi_u'),
                'meter second-1',
                'vertically integrated u-momentum component',
                time_name='uclm_time'
            )
    
            create_field(
                nc, 'vbar',
                ('vclm_time', 'eta_v', 'xi_v'),
                'meter second-1',
                'vertically integrated v-momentum component',
                time_name='vclm_time'
            )
    
            # Sea surface height 
            create_field(
                nc, 'SSH',
                ('ssh_time', 'eta_rho', 'xi_rho'),
                'meter',
                'sea surface height',
                time_name='ssh_time',
                coordinates='lon_rho lat_rho ssh_time'
            )
    
            create_field(
                nc, 'zeta',
                ('zeta_time', 'eta_rho', 'xi_rho'),
                'meter',
                'sea surface height',
                time_name='zeta_time',
                coordinates='lon_rho lat_rho zeta_time'
            )

    def create_grid_nc(self,output_file, inputs, outputs,prt_grd=None):    
        """
        Create and save a new CROCO grid file
        """
#        prt_grd=[AGRIF,prt_file,coef,imi,imax,jmin,jmax]
        if prt_grd is not None:
            if prt_grd[0]==True: # Means we are in AGRIF
                lev=prt_grd[1][-1]
                if not lev.isnumeric():
                    grid_name='croco_grd.nc.1'
                else:
                    grid_name=prt_grd[1][0:-1]+str(int(lev)+1)
                output_file=output_file+grid_name
		 
        nc = netcdf.Dataset(output_file, 'w', format='NETCDF4')

        # create global variables
        nc.created = datetime.now().isoformat()
        nc.type = 'CROCO grid file produced by easygrid_python.py'
        nc.createDimension('one', 1)

        if prt_grd is not None and prt_grd[0]: #AGRIF case
            nc.nx=outputs.h.shape[1]-2
            nc.ny=outputs.h.shape[0]-2

            # create dimensions
            nc.createDimension('xi_rho', outputs.h.shape[1])
            nc.createDimension('eta_rho', outputs.h.shape[0])
            nc.createDimension('xi_u', outputs.h.shape[1] - 1)
            nc.createDimension('eta_v', outputs.h.shape[0] - 1)
            nc.createDimension('xi_psi', outputs.h.shape[1] - 1)
            nc.createDimension('eta_psi', outputs.h.shape[0] - 1)
            nc.createDimension('four', 4)

            # Some empty variables in AGRIF
            nc.createVariable('xl', 'f8', ('one'))
            nc.variables['xl'].long_name = 'domain length in the XI-direction'
            nc.variables['xl'].units = 'meters'

            nc.createVariable('el', 'f8', ('one'))
            nc.variables['el'].long_name = 'domain length in the ETA-direction'
            nc.variables['el'].units = 'meters'
            nc.variables['el'][:] = inputs.ny

        else: # Usual case
            nc.nx = np.int32(inputs.nx)
            nc.ny = np.int32(inputs.ny)
            nc.size_x = inputs.size_x
            nc.size_y = inputs.size_y
            nc.tra_lon = inputs.tra_lon
            nc.tra_lat = inputs.tra_lat
            nc.rotation = inputs.rot
		    
            # create dimensions
            nc.createDimension('xi_rho', inputs.nx + 2)
            nc.createDimension('eta_rho', inputs.ny + 2)
            nc.createDimension('xi_u', inputs.nx + 1)
            nc.createDimension('eta_v', inputs.ny + 1)
            nc.createDimension('xi_psi', inputs.nx + 1)
            nc.createDimension('eta_psi', inputs.ny + 1)

            nc.createVariable('xl', 'f8', ('one'))
            nc.variables['xl'].long_name = 'domain length in the XI-direction'
            nc.variables['xl'].units = 'meters'
            nc.variables['xl'][:] = inputs.nx

            nc.createVariable('el', 'f8', ('one'))
            nc.variables['el'].long_name = 'domain length in the ETA-direction'
            nc.variables['el'].units = 'meters'
            nc.variables['el'][:] = inputs.ny

		
	# create variables and attributes
        nc.createVariable('spherical', 'S1', ('one'))
        nc.variables['spherical'].long_name = 'Grid type logical switch'
        nc.variables['spherical'].option_T = 'spherical'
        nc.variables['spherical'][:] = 'T'
		
        nc.createVariable('angle', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['angle'].long_name = 'angle between xi axis and east'
        nc.variables['angle'].units = 'radians' 
        nc.variables['angle'][:] = outputs.angle

        nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['h'].long_name = 'Final bathymetry at RHO-points'
        nc.variables['h'].units = 'meter'
        nc.variables['h'][:] = outputs.h

        nc.createVariable('hraw', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['hraw'].long_name = 'Working bathymetry at RHO-points'
        nc.variables['hraw'].units = 'meter'
        nc.variables['hraw'][:] = outputs.hraw

        nc.createVariable('f', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['f'].long_name = 'Coriolis parameter at RHO-points'
        nc.variables['f'].units = 'second-1'
        nc.variables['f'][:] = (4 * np.pi * np.sin(np.deg2rad(outputs.lat_rho)) /
		                        (23.9344699 * 3600))

        nc.createVariable('pm', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pm'].long_name = 'curvilinear coordinate metric in XI'
        nc.variables['pm'].units = 'meter-1'
        nc.variables['pm'][:] = outputs.pm

        nc.createVariable('pn', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pn'].long_name = 'curvilinear coordinate metric in ETA'
        nc.variables['pn'].units = 'meter-1'
        nc.variables['pn'][:] = outputs.pn

        nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
        nc.variables['lon_rho'].units = 'degree_east'
        nc.variables['lon_rho'][:] = outputs.lon_rho

        nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
        nc.variables['lat_rho'].units = 'degree_north'
        nc.variables['lat_rho'][:] = outputs.lat_rho

        nc.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['mask_rho'].long_name = 'mask on RHO-points'
        nc.variables['mask_rho'].option_0 = 'land'
        nc.variables['mask_rho'].option_1 = 'water'
        nc.variables['mask_rho'][:] = outputs.mask_rho

    	# Extraneous variables should be placed at the end (ensures no
    	# later problems with e.g., partit
        nc.createVariable('lon_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lon_psi'].long_name = 'longitude of PSI-points'
        nc.variables['lon_psi'].units = 'degree_east'
        if prt_grd is not None and prt_grd[0]:
            nc.variables['lon_psi'][:] = outputs.lon_psi
        else:
            nc.variables['lon_psi'][:] = outputs.lon_psi[1:-1, 1:-1]

        nc.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lat_psi'].long_name = 'latitude of PSI-points'
        nc.variables['lat_psi'].units = 'degree_north'
        if prt_grd is not None and prt_grd[0]:
            nc.variables['lat_psi'][:] = outputs.lat_psi
        else:
            nc.variables['lat_psi'][:] = outputs.lat_psi[1:-1, 1:-1]

        nc.createVariable('lon_u', 'f8', ('eta_rho', 'xi_u'))
        nc.variables['lon_u'].long_name = 'longitude of U-points'
        nc.variables['lon_u'].units = 'degree_east'
        nc.variables['lon_u'][:] = outputs.lon_u

        nc.createVariable('lat_u', 'f8', ('eta_rho', 'xi_u'))
        nc.variables['lat_u'].long_name = 'latitude of U-points'
        nc.variables['lat_u'].units = 'degree_north'
        nc.variables['lat_u'][:] = outputs.lat_u

        nc.createVariable('lon_v', 'f8', ('eta_v', 'xi_rho'))
        nc.variables['lon_v'].long_name = 'longitude of V-points'
        nc.variables['lon_v'].units = 'degree_east'
        nc.variables['lon_v'][:] = outputs.lon_v

        nc.createVariable('lat_v', 'f8', ('eta_v', 'xi_rho'))
        nc.variables['lat_v'].long_name = 'latitude of RHO-points'
        nc.variables['lat_v'].units = 'degree_north'
        nc.variables['lat_v'][:] = outputs.lat_v

        if prt_grd is not None and prt_grd[0]: 
            nc.createVariable('refine_coef', 'i', ('one'))
            nc.variables['refine_coef'].long_name ='Grid refinement coefficient'
            nc.variables['refine_coef'][:]=prt_grd[2]

            nc.createVariable('grd_pos','i',('four'))
            nc.variables['grd_pos'].long_name='Subgrid location in the parent grid: psi corner points (imin imax jmin jmax)'
            nc.variables['grd_pos'][:]=prt_grd[3:]


        nc.close()
        print('Writting '+output_file+' done')

        if prt_grd is not None and prt_grd[0]:
            print('Create an AGRIF_FixedGrids.in file')
            fname='AGRIF_FixedGrids.in'
            fid=open(fname,'w')
            fid.write('    1\n')#'%s\n','    1');
            fid.write('    '+str(prt_grd[3])+ \
                               '    '+str(prt_grd[4])+ \
                               '    '+str(prt_grd[5])+\
                               '    '+str(prt_grd[6])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2])+\
                               '    '+str(prt_grd[2]))
            fid.write('\n    0')
            fid.write('\n# number of children per parent')
            fid.write('\n# imin imax jmin jmax spacerefx spacerefy timerefx timerefy')
            fid.write('\n# [all coordinates are relative to each parent grid!]')
            fid.write('\n~')
            fid.close()

    def create_tide_nc(self,filename, grdobj, created_by='make_tides.py',cur=False,pot=False):
        # Global attributes
        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type  = 'CROCO tide file produced by %s' % created_by

        # Dimensions
        nc.createDimension('xi_rho', grdobj.lon.shape[1])
        nc.createDimension('eta_rho', grdobj.lon.shape[0])
        nc.createDimension('tide_period', None)
        
        # Create the variables and write...
        nc.createVariable('tide_period', 'f8', ('tide_period'))
        nc.variables['tide_period'].long_name = 'Tide angular period'
        nc.variables['tide_period'].units = 'Hours'

        nc.createVariable('tide_Ephase', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
        nc.variables['tide_Ephase'].long_name = 'Tidal elevation phase angle'
        nc.variables['tide_Ephase'].units = 'Degrees' 

        nc.createVariable('tide_Eamp', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
        nc.variables['tide_Eamp'].long_name = 'Tidal elevation amplitude'
        nc.variables['tide_Eamp'].units = 'Meter'

        if cur:

            nc.createVariable('tide_Cmin', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Cmin'].long_name = 'Tidal current ellipse semi-minor axis'
            nc.variables['tide_Cmin'].units = 'm.s-1'

            nc.createVariable('tide_Cmax', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Cmax'].long_name = 'Tidal current, ellipse semi-major axis'
            nc.variables['tide_Cmax'].units = 'm.s-1'

            nc.createVariable('tide_Cangle', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Cangle'].long_name = 'Tidal current inclination angle'
            nc.variables['tide_Cangle'].units = 'Degrees'

            nc.createVariable('tide_Cphase', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Cphase'].long_name = 'Tidal current phase angle'
            nc.variables['tide_Cphase'].units = 'Meter'

        if pot:

            nc.createVariable('tide_Pamp', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Pamp'].long_name = 'Tidal potential amplitude'
            nc.variables['tide_Pamp'].units = 'Meter'

            nc.createVariable('tide_Pphase', 'f8',('tide_period', 'eta_rho', 'xi_rho'))
            nc.variables['tide_Pphase'].long_name = 'Tidal potential phase angle'
            nc.variables['tide_Pphase'].units = 'Degrees'

        nc.close()

    def create_river_nc(self,filename, grdobj,Nsrc,TS,created_by='make_river.py'):

        nc = netcdf.Dataset(filename, 'w', format='NETCDF4')
        nc.created = datetime.now().isoformat()
        nc.type  = 'CROCO river file produced by %s' % created_by
        nc.grd_file = grdobj.grid_file
        
        # Dimensions
        nc.createDimension('qbar_time',0)
        nc.createDimension('n_qbar',Nsrc)

        # Create the variables and write...
        nc.createVariable('runoff_name','S30', ('n_qbar',))
        nc.variables['runoff_name'].long_name = 'River Name'

        nc.createVariable('qbar_time',np.float64, ('qbar_time',))
        nc.variables['qbar_time'].long_name = 'river discharge time'

        nc.createVariable('Qbar',np.float64, ('n_qbar','qbar_time',))
        nc['Qbar'].long_name = 'river discharge'
        nc['Qbar'].units = 'm3 s-1'
   
        if TS:
            nc.createDimension('temp_time',0)
            nc.createDimension('salt_time',0)

            nc.createVariable('temp_src_time',np.float64, ('temp_time',))
            nc.variables['temp_src_time'].long_name = 'river temperature time'

            nc.createVariable('salt_src_time',np.float64, ('salt_time',))
            nc.variables['salt_src_time'].long_name = 'river salinity time'
    
            nc.createVariable('temp_src',np.float64, ('n_qbar','temp_time',))
            nc['temp_src'].long_name = 'river temperature'
            nc['temp_src'].units = 'Celsius'

            nc.createVariable('salt_src',np.float64, ('n_qbar','salt_time',))
            nc['salt_src'].long_name = 'river salinity'
            nc['salt_src'].units = '1'

        nc.close()











