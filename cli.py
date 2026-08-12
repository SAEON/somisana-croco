'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow python functions to be run from the cli docker image for this repo
So this is the entry point for the docker image (see Dockerfile.cli).
But it's also handy if you want to execute a python function from inside a bash script
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import argparse
import sys, os
from datetime import datetime, timedelta
import calendar
from crocotools_py.preprocess import make_tides,reformat_gfs_atm,reformat_saws_atm,make_ini,make_bry,make_clm
from crocotools_py.postprocess import (get_ts_multivar, croco_srf_2_ww3,
                                        read_timeseries_coords)
from crocotools_py.plotting import plot as crocplot
from crocotools_py.regridding import regrid_tier1, regrid_tier2, regrid_tier3
from crocotools_py.products import (make_products_base, add_anomalies, add_mhw_mcs,
                                    add_sst_front, add_stratification)
from crocotools_py.plotting_products import (plot_coastal_flags, plot_timeseries_mhw,
                                             plot_timeseries_stratification)

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYY-MM-DD HH:MM:SS'.")

def parse_int(value):
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: {value}")

def parse_list(value):
    return [float(x) for x in value.split(',')]

def parse_list_str(value):
    if value is None or value == 'None':
        return None
    else:
        return [x.strip() for x in value.split(',')]
    
def parse_bool(s: str) -> bool:
    try:
        return {'true':True, 'false':False}[s.lower()]
    except KeyError:
        raise argparse.ArgumentTypeError(f'expect true/false, got: {s}')

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-croco repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...

    # ------------------
    # reformat_saws_atm
    # ------------------
    parser_reformat_saws_atm = subparsers.add_parser('reformat_saws_atm', 
            help='convert the SAWS UM nc files into nc files which can be ingested by CROCO using the ONLINE cpp key')
    parser_reformat_saws_atm.add_argument('--sawsDir', required=True, help='Directory containing the SAWS UM files')
    parser_reformat_saws_atm.add_argument('--backupDir', required=True, help='Directory containing already reformatted data, used for variables not provided by SAWS')
    parser_reformat_saws_atm.add_argument('--outputDir', required=True, help='Directory to save reformated nc files')
    parser_reformat_saws_atm.add_argument('--run_date', required=False, type=parse_datetime,
            default=None,
            help='initialisation time in format "YYYY-MM-DD HH:MM:SS"')
    parser_reformat_saws_atm.add_argument('--hdays', required=False, type=float, 
            default=5.,
            help='hindcast days i.e before run_date')
    parser_reformat_saws_atm.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model - reformatted file time will be in days since 1-Jan-Yorig')
    def reformat_saws_atm_handler(args):
        reformat_saws_atm(args.sawsDir,args.backupDir,args.outputDir,args.run_date,args.hdays,args.Yorig)
    parser_reformat_saws_atm.set_defaults(func=reformat_saws_atm_handler)
     
    # ------------------
    # reformat_gfs_atm
    # ------------------
    parser_reformat_gfs_atm = subparsers.add_parser('reformat_gfs_atm', 
            help='convert the downloaded GFS grb files into nc files which can be ingested by CROCO using the ONLINE cpp key')
    parser_reformat_gfs_atm.add_argument('--gfsDir', required=True, help='Directory containing the grb files downloaded using download_gfs_atm')
    parser_reformat_gfs_atm.add_argument('--outputDir', required=True, help='Directory to save reformated nc files')
    parser_reformat_gfs_atm.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model - reformatted file time will be in days since 1-Jan-Yorig')
    def reformat_gfs_atm_handler(args):
        reformat_gfs_atm(args.gfsDir,args.outputDir,args.Yorig)
    parser_reformat_gfs_atm.set_defaults(func=reformat_gfs_atm_handler)
    
    # --------------------------
    # plot/animate croco output
    # --------------------------
    parser_crocplot = subparsers.add_parser('crocplot', 
            help='do a plot/animation of a croco output file(s)')
    parser_crocplot.add_argument('--fname', required=True, type=str, help='input native CROCO filename (can handle wildcards to animate over multiple files)')
    parser_crocplot.add_argument('--var', required=False, default='temp', type=str, help='the variable name to plot')
    parser_crocplot.add_argument('--gif_out', required=False, type=str, help='the output gif filename')
    parser_crocplot.add_argument('--mp4_out', required=False, type=str, help='the output mp4 filename')
    parser_crocplot.add_argument('--level', required=False, default=None, type=parse_int, help='level to plot. If >=0, then a sigma level is plotted. If <0 then a z level (in m) is plotted. Default behaviour will plot the surface layer')
    parser_crocplot.add_argument('--ticks', required=False, type=parse_list,
                         default=None, 
                         help='contour ticks to use in plotting the variable')
    parser_crocplot.add_argument('--cbar_label', required=False, default=r'temperature ($\degree$C)', type=str, help='the label used for the colorbar')
    parser_crocplot.add_argument('--cmap', required=False, default='Spectral_r', type=str,
            help="the matplotlib colormap to use. 'mhw_mcs' is registered by crocotools_py.plotting for marine heatwave / cold spell categories - use it with --ticks -4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5 to get one colour band per category")
    parser_crocplot.add_argument('--isobaths', required=False, type=parse_list,
                         default=[100,500],
                         help='the isobaths to add to the figure')
    parser_crocplot.add_argument('--add_vectors', type=parse_bool, 
                       default='True',
                       help='If True, current vectors are added to the plot')
    parser_crocplot.add_argument('--skip_time', required=False, type=int, default=1,
            help='Number of time-steps to skip between frames in the animation')
    parser_crocplot.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    def crocplot_handler(args):
        # a lot of the crocplot inputs are hard coded below but could be made configurable in future
        # there are also other potential optional inputs to this function which we aren't specifying here
        # so this is a work in progress...
        crocplot(args.fname, var=args.var,
                      grdname=None, # could make this configurable in the cli args
                      time=slice(None),
                      level=args.level,
                      ticks = args.ticks,
                      cmap = args.cmap,
                      extents = None,
                      Yorig = args.Yorig,
                      cbar_label=args.cbar_label,
                      add_vectors = args.add_vectors,
                      skip_time = args.skip_time,
                      isobaths=args.isobaths,
                      gif_out=args.gif_out,
                      mp4_out=args.mp4_out
                      )
    parser_crocplot.set_defaults(func=crocplot_handler)
    
    # --------------
    # regrid_tier1
    # --------------
    parser_regrid_tier1 = subparsers.add_parser('regrid_tier1', 
            help='tier 1 regridding of a raw CROCO output file: regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid -> rotates u/v to be east/north components instead of grid-aligned components -> adds a depth variable providing the depths of each sigma level at each time-step')
    parser_regrid_tier1.add_argument('--fname', required=True, type=str, help='input native CROCO filename - can include wildcards (*) to process multiple files in a single command')
    parser_regrid_tier1.add_argument('--grdname', required=False, type=str,
                                     help='optional input to provide a separate grid file (if grid vars arent in your output files)')
    parser_regrid_tier1.add_argument('--dir_out', required=True, help='tier 1 output directory')
    parser_regrid_tier1.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_regrid_tier1.add_argument('--doi_link', required=False, type=str, help='Doi link to where the data can be located.')
    def regrid_tier1_handler(args):
        regrid_tier1(args.fname, args.dir_out, args.grdname, args.Yorig, args.doi_link)
    parser_regrid_tier1.set_defaults(func=regrid_tier1_handler)
    
    # --------------
    # regrid_tier2
    # --------------
    parser_regrid_tier2 = subparsers.add_parser('regrid_tier2', 
            help='tier 2 regridding of a raw CROCO output file: regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid -> rotates u/v to be east/north components instead of grid-aligned components -> regrids the sigma levels to the user defined constant z levels, including the surface and bottom layers -> output variables are the same as tier 1, only depths is now a dimension with the user specified values')
    parser_regrid_tier2.add_argument('--fname', required=True, type=str, help='input native CROCO filename - can include wildcards (*) to process multiple files in a single command')
    parser_regrid_tier2.add_argument('--grdname', required=False, type=str,
                                     help='optional input to provide a separate grid file (if grid vars arent in your output files)')
    parser_regrid_tier2.add_argument('--dir_out', required=True, help='tier 2 output directory')
    parser_regrid_tier2.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_regrid_tier2.add_argument('--doi_link', required=False, type=str, help='Doi link to where the data can be located.')
    parser_regrid_tier2.add_argument('--depths', required=False, type=parse_list,
                         default=[0,-5,-10,-20,-50,-100,-200,-500,-1000],
                         help='list of depths to extract (in metres, negative down)')
    parser_regrid_tier2.add_argument('--varList', required=False, type=parse_list_str,
                         default=['temp','salt','u','v'],
                         help='comma separated list of variables to regrid e.g. temp,salt,u,v')
    def regrid_tier2_handler(args):
        regrid_tier2(args.fname, args.dir_out, args.grdname, args.Yorig, args.doi_link, depths = args.depths, varList = args.varList)
    parser_regrid_tier2.set_defaults(func=regrid_tier2_handler)
    
    # --------------
    # regrid_tier3
    # --------------
    parser_regrid_tier3 = subparsers.add_parser('regrid_tier3', 
            help='tier 3 regridding of a CROCO output: takes the output of regrid-tier2 as input and regrids the horizontal grid to a regular grid with a specified grid spacing. Output variables are the same as tier 1 and 2, only horizontal grid is now rectilinear with hz dimensions of longitude,latitude i.e. horizontal grid is no longer curvilinear. The extents of the rectilinear grid are automatically determined using the curvilinear grid extents.')
    parser_regrid_tier3.add_argument('--fname', required=True, type=str, help='input regridded tier2 filename - can include wildcards (*) to process multiple files in a single command')
    parser_regrid_tier3.add_argument('--dir_out', required=True, help='tier 3 output directory')
    parser_regrid_tier3.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_regrid_tier3.add_argument('--doi_link', required=False, type=str, help='Doi link to where the data can be located.')
    parser_regrid_tier3.add_argument('--spacing', type=float,default=0.01,
                         help='constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output')
    parser_regrid_tier3.add_argument('--method', type=str, default='nearest',
                         choices=['nearest', 'linear', 'cubic'],
                         help="interpolation method passed to scipy.interpolate.griddata (default='nearest')")
    def regrid_tier3_handler(args):
        regrid_tier3(args.fname, args.dir_out, args.Yorig, args.doi_link,
                     spacing=args.spacing, method=args.method)
    parser_regrid_tier3.set_defaults(func=regrid_tier3_handler)
    
    # ----------------
    # get_ts_multivar
    # ----------------
    parser_get_ts_multivar = subparsers.add_parser('get_ts_multivar',
            help='extract time-series (or profiles through time) of several variables at one or many sites, into a single netcdf file with a site dimension')
    parser_get_ts_multivar.add_argument('--fname', required=True, type=str,
            help='input CROCO filename - can include wildcards (*) to extract across multiple files')
    parser_get_ts_multivar.add_argument('--coords_file', required=False, type=str, default=None,
            help="csv text file of the sites to extract, one per line as 'name,lon,lat' (see configs/<domain>/timeseries_coords.csv). Use this, or --lon/--lat for a single unnamed site")
    parser_get_ts_multivar.add_argument('--lon', required=False, type=float, default=None,
            help='longitude of a single extraction site, if --coords_file is not used')
    parser_get_ts_multivar.add_argument('--lat', required=False, type=float, default=None,
            help='latitude of a single extraction site, if --coords_file is not used')
    parser_get_ts_multivar.add_argument('--varList', required=False, type=parse_list_str,
                        default=['temp','salt','u','v'],
                        help="comma separated list of CROCO variable names e.g. 'temp,salt,u,v'. u/v type variables can only be rotated to east/north components as a pair, so naming either component extracts both")
    parser_get_ts_multivar.add_argument('--Yorig', type=parse_int,
                        default=2000,
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be in seconds since Yorig-01-01')
    parser_get_ts_multivar.add_argument('--fname_out', required=True, help='output filename')
    def get_ts_multivar_handler(args):
        if args.coords_file is not None:
            site_names, lon, lat = read_timeseries_coords(args.coords_file)
        elif args.lon is not None and args.lat is not None:
            site_names, lon, lat = None, args.lon, args.lat
        else:
            raise SystemExit('get_ts_multivar needs either --coords_file, or both --lon and --lat')
        get_ts_multivar(args.fname, lon, lat,
               Yorig=args.Yorig,
               varList=args.varList,
               site_names=site_names,
               nc_out=args.fname_out)
    parser_get_ts_multivar.set_defaults(func=get_ts_multivar_handler)
    
    # ----------------
    # make_tides_fcst
    # ----------------
    parser_make_tides_fcst = subparsers.add_parser('make_tides_fcst',
            help='Make a tidal forcing file for as part of the CROCO operational workflow')
    parser_make_tides_fcst.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the raw tidal files from e.g. TPXO')
    parser_make_tides_fcst.add_argument('--output_dir', required=True, type=str,
            help='Path to where the forcing file will be saved. This directory also needs a crocotools_param.py file')
    parser_make_tides_fcst.add_argument('--run_date', required=True, type=parse_datetime, 
            help='operational workflow initialisation time in format "YYYY-MM-DD HH:MM:SS"')
    parser_make_tides_fcst.add_argument('--hdays', required=True, type=int,
            help='Number of days before run_date, corresponding to the run start time')
    parser_make_tides_fcst.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model')
    def make_tides_fcst_handler(args):
        
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        fname_out = params.croco_prefix + args.run_date.strftime('_%Y%m%d_%H.nc')
        run_date_ini = args.run_date - timedelta(days=args.hdays)
        
        make_tides(args.input_dir,args.output_dir,run_date_ini,args.Yorig,fname_out)
        
    parser_make_tides_fcst.set_defaults(func=make_tides_fcst_handler)
    
    # ----------------
    # make_tides_inter
    # ----------------
    parser_make_tides_inter = subparsers.add_parser('make_tides_inter',
            help='Make monthly tidal forcing files for CROCO interannual runs')
    parser_make_tides_inter.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the raw tidal files from e.g. TPXO')
    parser_make_tides_inter.add_argument('--output_dir', required=True, type=str,
            help='Path to where the forcing file will be saved. This directory also needs a crocotools_param.py file')
    parser_make_tides_inter.add_argument('--month_start', required=True, type=str, 
            help='first month in the interannual run in format "YYYY-MM"')
    parser_make_tides_inter.add_argument('--month_end', required=True, type=str,
            help='last month in the interannual run in format "YYYY-MM"')
    parser_make_tides_inter.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model')
    def make_tides_inter_handler(args):
        
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        month_now = datetime.strptime(args.month_start+'-01','%Y-%m-%d')
        month_end = datetime.strptime(args.month_end+'-01','%Y-%m-%d')
        
        while month_now <= month_end:
            
            print('working on '+month_now.strftime('%Y-%m'))
            
            fname_out = params.croco_prefix + month_now.strftime('_Y%YM%m.nc') + params.croco_suffix
            
            make_tides(args.input_dir,args.output_dir,month_now,args.Yorig,fname_out)
            
            month_now=month_now+timedelta(days=32) # 32 days ensures we get to the next month
            month_now=datetime(month_now.year, month_now.month, 1) # set to the first day of the month
        
    parser_make_tides_inter.set_defaults(func=make_tides_inter_handler)
    
    # ----------------
    # make_ini_fcst
    # ----------------
    parser_make_ini_fcst = subparsers.add_parser('make_ini_fcst',
            help='Make initial condition file from OGCM data as part of the CROCO operational workflow.')
    parser_make_ini_fcst.add_argument('--input_file', required=True, type=str, 
            help='Path and filename of the OGCM input file i.e. "path/to/file/direcory/and/filename.nc"')
    parser_make_ini_fcst.add_argument('--output_dir', required=True, type=str,
            help='Path to where the ini file will be saved. This directory also needs a crocotools_param.py file')
    parser_make_ini_fcst.add_argument('--run_date', required=True, type=parse_datetime, 
            help='operational workflow initialisation time in format "YYYY-MM-DD HH:MM:SS"')
    parser_make_ini_fcst.add_argument('--hdays', required=True, type=int,
            help='Number of days before run_date to initialise model')
    parser_make_ini_fcst.add_argument('--Yorig', required=True, type=int,
            help='the Yorig value used in setting up the CROCO model')
    def make_ini_fcst_handler(args):
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        fname_out = params.ini_prefix + args.run_date.strftime('_%Y%m%d_%H.nc')
        ini_date = args.run_date - timedelta(days=args.hdays)
        
        make_ini(args.input_file,args.output_dir,ini_date,args.Yorig,fname_out)
    
    parser_make_ini_fcst.set_defaults(func=make_ini_fcst_handler)
    
    # ----------------
    # make_ini_inter
    # ----------------
    parser_make_ini_inter = subparsers.add_parser('make_ini_inter',
            help='Make initial condition file from OGCM data for an inter-annual run.')
    parser_make_ini_inter.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the monthly OGCM files')
    parser_make_ini_inter.add_argument('--output_dir', required=True, type=str,
            help='Path to where the ini file will be saved. This directory also needs a crocotools_param.py file')
    parser_make_ini_inter.add_argument('--month_start', required=True, type=str, 
            help='first month in the interannual run in format "YYYY-MM"')
    parser_make_ini_inter.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model')
    def make_ini_inter_handler(args):
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        ini_date = datetime.strptime(args.month_start+'-01','%Y-%m-%d')
        fname_in = os.path.join(args.input_dir, ini_date.strftime(params.input_file_fmt))
        
        fname_out = params.ini_prefix + ini_date.strftime('_Y%YM%m.nc') + params.ini_suffix
        
        make_ini(fname_in,args.output_dir,ini_date,args.Yorig,fname_out)
    
    parser_make_ini_inter.set_defaults(func=make_ini_inter_handler)
    
    # ----------------
    # make_bry_fcst
    # ----------------
    parser_make_bry_fcst = subparsers.add_parser('make_bry_fcst',
            help='Make ocean boundary conditions as part of the CROCO operational workflow.')
    parser_make_bry_fcst.add_argument('--input_file', required=True, type=str, 
            help='Path and filename of input file i.e. "path/to/file/direcory/and/filename.nc"')
    parser_make_bry_fcst.add_argument('--output_dir', required=True, type=str,
            help='Path to where the bry file will be saved. This directory also needs a crocotools_param.py file')
    parser_make_bry_fcst.add_argument('--run_date', required=True, type=parse_datetime, 
            help='Operational workflow initialisation time in format "YYYY-MM-DD HH:MM:SS"')
    parser_make_bry_fcst.add_argument('--hdays', required=True, type=int,
            help='Number of hindcast days to add to the boundary file')
    parser_make_bry_fcst.add_argument('--fdays', required=True, type=int,
            help='Number of forecast days to add to the boundary file')
    parser_make_bry_fcst.add_argument('--Yorig', required=True, type=int,
                        help='the Yorig value used in setting up the CROCO model')
    def make_bry_fcst_handler(args):
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        # create a 2 day buffer around the time span of the CROCO simulation
        # (2 days is just to be safe - the nearest available times to ini_date and end_date are used in make_bry())
        hdays = args.hdays + 2
        fdays = args.fdays + 2
        
        fname_out = params.bry_prefix + args.run_date.strftime('_%Y%m%d_%H.nc')
        ini_date = args.run_date - timedelta(days=hdays)
        end_date = args.run_date + timedelta(days=fdays)
        
        make_bry(args.input_file,args.output_dir,ini_date,end_date,args.Yorig,fname_out)
        
    parser_make_bry_fcst.set_defaults(func=make_bry_fcst_handler)
    
    # ----------------
    # make_bry_inter
    # ----------------
    parser_make_bry_inter = subparsers.add_parser('make_bry_inter',
            help='Make monthly ocean boundary condition files for CROCO interannual runs')
    parser_make_bry_inter.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the monthly OGCM files')
    parser_make_bry_inter.add_argument('--output_dir', required=True, type=str,
            help='Path to where the forcing files will be saved. This directory also needs a crocotools_param.py file')
    parser_make_bry_inter.add_argument('--month_start', required=True, type=str, 
            help='first month in the interannual run in format "YYYY-MM"')
    parser_make_bry_inter.add_argument('--month_end', required=True, type=str,
            help='last month in the interannual run in format "YYYY-MM"')
    parser_make_bry_inter.add_argument('--Yorig', required=True, type=int,
            help='the Yorig value used in setting up the CROCO model')
    def make_bry_inter_handler(args):
        
        sys.path.append(args.output_dir)
        import crocotools_param as params
        
        month_now = datetime.strptime(args.month_start+'-01','%Y-%m-%d')
        month_end = datetime.strptime(args.month_end+'-01','%Y-%m-%d')
        
        while month_now <= month_end:
            
            print('working on '+month_now.strftime('%Y-%m'))
            
            # define a list of 3 input files - last month, this month and next month
            # This is to ensure that we always have boundary data for the start and end of the CROCO run for this month
            month_prev = month_now - timedelta(days=10) # an arbitrary date in the pervious month
            month_next = month_now + timedelta(days=32) # 32 days ensures we get to the next month
            fname_month_prev = os.path.join(args.input_dir, month_prev.strftime(params.input_file_fmt))
            fname_month_now = os.path.join(args.input_dir, month_now.strftime(params.input_file_fmt))
            fname_month_next = os.path.join(args.input_dir, month_next.strftime(params.input_file_fmt))
            # We could handle the case where the user doesn't have files for the previous and next month
            # (to do this I think we'd need to create the input_file list with file names which do exist,
            # and then we'd need a check inside make_bry() where we check that ini_date and end_date are 
            # covered by the OGCM files, and if not we'd pad with the nearest available values)
            # But I'm in a hurry so for now we'll just make sure that these files do exist
            if not os.path.exists(fname_month_prev):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_prev)
            if not os.path.exists(fname_month_now):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_now)
            if not os.path.exists(fname_month_next):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_next)
            input_file=[fname_month_prev,
                        fname_month_now,
                        fname_month_next]
            
            # define ini_date and end_date using a 1 day buffer either side of the month
            ini_date = month_now - timedelta(days=1)
            day_end = calendar.monthrange(month_now.year,month_now.month)[1]
            end_date = datetime(month_now.year,month_now.month,day_end) + timedelta(days=2)    
            
            # make the boundary file for this month
            fname_out = params.bry_prefix + month_now.strftime('_Y%YM%m.nc')
            make_bry(input_file,args.output_dir,ini_date,end_date,args.Yorig,fname_out)
            
            month_now=datetime(month_next.year, month_next.month, 1) # set month_now to the first day of the next month
        
    parser_make_bry_inter.set_defaults(func=make_bry_inter_handler)

    # ----------------------
    # make_products_base
    # ----------------------
    parser_products_base = subparsers.add_parser('make_products_base',
            help='Create the daily averaged products file from raw CROCO output. The output is itself a valid CROCO file (it carries the full grid), so everything in postprocess.py works on it. This is the file which the add_* commands then append their products to')
    parser_products_base.add_argument('--fname_in', required=True, type=str,
            help='input raw CROCO filename (or a file pattern to be used with open_mfdataset())')
    parser_products_base.add_argument('--fname_out', required=True, type=str,
            help='output filename, conventionally croco_avg_products.nc')
    parser_products_base.add_argument('--Yorig', required=False, type=parse_int, default=2000,
            help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_products_base.add_argument('--varList', required=False, type=parse_list_str, default=None,
            help="comma separated list of variables to daily-average, e.g. 'temp,salt' (the default). 'zeta' is always included, as it is needed to compute the depths of the sigma levels")
    def make_products_base_handler(args):
        make_products_base(args.fname_in, args.fname_out,
                           Yorig=args.Yorig,
                           varList=args.varList)
    parser_products_base.set_defaults(func=make_products_base_handler)

    # ----------------------
    # add_anomalies
    # ----------------------
    parser_add_anom = subparsers.add_parser('add_anomalies',
            help='Add daily temperature and sea level anomalies to the products file, computed against the day-of-year climatology (the same baseline used for the MHW/MCS categories)')
    parser_add_anom.add_argument('--fname_out', required=True, type=str,
            help='the products file to add to, as created by make_products_base. If it does not exist it gets created first, in which case --fname is needed')
    parser_add_anom.add_argument('--clim_file', required=True, type=str,
            help='pre-built day-of-year climatology file (see the hindcasts/*/climatology directory for the domain)')
    parser_add_anom.add_argument('--fname_in', required=False, type=str, default=None,
            help='CROCO file(s) to read the daily fields from - can include wildcards (*). Defaults to fname_out, which is the usual case of adding to a products file that already exists. If fname_out does not exist yet it gets created from this input first')
    parser_add_anom.add_argument('--Yorig', required=False, type=parse_int, default=2000,
            help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    def add_anomalies_handler(args):
        add_anomalies(args.fname_out, args.clim_file, fname_in=args.fname_in, Yorig=args.Yorig)
    parser_add_anom.set_defaults(func=add_anomalies_handler)

    # ----------------------
    # add_mhw_mcs
    # ----------------------
    parser_add_mhw = subparsers.add_parser('add_mhw_mcs',
            help='Add Marine Heatwave (MHW) and Marine Cold Spell (MCS) event categories to the products file, using a pre-built day-of-year climatology and percentile thresholds. Categories are signed: +1..+4 = MHW, 0 = no event, -1..-4 = MCS')
    parser_add_mhw.add_argument('--fname_out', required=True, type=str,
            help='the products file to add to, as created by make_products_base. If it does not exist it gets created first, in which case --fname is needed')
    parser_add_mhw.add_argument('--clim_file', required=True, type=str,
            help='pre-built day-of-year climatology file (see the hindcasts/*/climatology directory for the domain)')
    parser_add_mhw.add_argument('--thresh_file', required=True, type=str,
            help='pre-built day-of-year percentile threshold file, containing threshold_90 and threshold_10')
    parser_add_mhw.add_argument('--fname_in', required=False, type=str, default=None,
            help='CROCO file(s) to read the daily fields from - can include wildcards (*). Defaults to fname_out, which is the usual case of adding to a products file that already exists. If fname_out does not exist yet it gets created from this input first')
    parser_add_mhw.add_argument('--Yorig', required=False, type=parse_int, default=2000,
            help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_add_mhw.add_argument('--batch_size', required=False, type=parse_int, default=5,
            help='number of eta_rho rows processed at a time (trade-off between memory use and speed)')
    def add_mhw_mcs_handler(args):
        add_mhw_mcs(args.fname_out, args.clim_file, args.thresh_file,
                    fname_in=args.fname_in, Yorig=args.Yorig, batch_size=args.batch_size)
    parser_add_mhw.set_defaults(func=add_mhw_mcs_handler)

    # ----------------------
    # add_sst_front
    # ----------------------
    parser_add_front = subparsers.add_parser('add_sst_front',
            help='Add the daily surface thermal front magnitude (the horizontal gradient of daily mean SST, in degC/km) to the products file')
    parser_add_front.add_argument('--fname_out', required=True, type=str,
            help='the products file to add to, as created by make_products_base. If it does not exist it gets created first, in which case --fname is needed')
    parser_add_front.add_argument('--fname_in', required=False, type=str, default=None,
            help='CROCO file(s) to read the daily fields from - can include wildcards (*). Defaults to fname_out, which is the usual case of adding to a products file that already exists. If fname_out does not exist yet it gets created from this input first')
    parser_add_front.add_argument('--Yorig', required=False, type=parse_int, default=2000,
            help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    def add_sst_front_handler(args):
        add_sst_front(args.fname_out, fname_in=args.fname_in, Yorig=args.Yorig)
    parser_add_front.set_defaults(func=add_sst_front_handler)

    # ----------------------
    # add_stratification
    # ----------------------
    parser_add_strat = subparsers.add_parser('add_stratification',
            help='Add the daily bottom density, the density at target_depth, and the difference between them (the stratification) to the products file')
    parser_add_strat.add_argument('--fname_out', required=True, type=str,
            help='the products file to add to, as created by make_products_base. If it does not exist it gets created first, in which case --fname is needed')
    parser_add_strat.add_argument('--fname_in', required=False, type=str, default=None,
            help='CROCO file(s) to read the daily fields from - can include wildcards (*). Defaults to fname_out, which is the usual case of adding to a products file that already exists. If fname_out does not exist yet it gets created from this input first')
    parser_add_strat.add_argument('--Yorig', required=False, type=parse_int, default=2000,
            help='Origin year used in setting up CROCO time i.e. CROCO time will be seconds since Yorig-01-01')
    parser_add_strat.add_argument('--target_depth', required=False, type=float, default=5.0,
            help='depth (in metres, positive down) whose density is compared against the bottom density to get the stratification')
    def add_stratification_handler(args):
        add_stratification(args.fname_out, fname_in=args.fname_in, Yorig=args.Yorig,
                           target_depth=args.target_depth)
    parser_add_strat.set_defaults(func=add_stratification_handler)

    # ----------------------
    # plot_coastal_flags
    # ----------------------
    parser_flags = subparsers.add_parser('plot_coastal_flags',
            help='Map the worst MHW/MCS category reached at each site over the forecast window, as coloured flags along the coast')
    parser_flags.add_argument('--ts_file', required=True, type=str,
            help='time-series file written by get_ts_multivar, containing the category variable')
    parser_flags.add_argument('--out_path', required=True, type=str,
            help='output png filename')
    parser_flags.add_argument('--level', required=False, type=parse_int, default=-1,
            help='sigma level to summarise: -1 (default) is the surface, 0 is the bottom')
    parser_flags.add_argument('--depth_name', required=False, type=str, default='Surface',
            help="label for the level, used in the title e.g. 'Surface' or 'Bottom'")
    parser_flags.add_argument('--title', required=False, type=str, default=None,
            help="first line of the plot title e.g. 'SA West Coast  ·  MHW / MCS Flag Map'")
    parser_flags.add_argument('--extent', required=False, type=parse_list, default=None,
            help='map extent as lon0,lon1,lat0,lat1. Derived from the site coordinates if not given')
    def plot_coastal_flags_handler(args):
        plot_coastal_flags(args.ts_file, args.out_path,
                           level=args.level,
                           depth_name=args.depth_name,
                           title=args.title,
                           extent=args.extent)
    parser_flags.set_defaults(func=plot_coastal_flags_handler)

    # ----------------------
    # plot_timeseries_mhw
    # ----------------------
    parser_ts_mhw = subparsers.add_parser('plot_timeseries_mhw',
            help='Plot per-site temperature time-series against the day-of-year climatology and the MHW/MCS thresholds')
    parser_ts_mhw.add_argument('--ts_file', required=True, type=str,
            help='time-series file written by get_ts_multivar, containing the temp variable')
    parser_ts_mhw.add_argument('--clim_file', required=True, type=str,
            help='pre-built day-of-year climatology file (as used as input to add_mhw_mcs)')
    parser_ts_mhw.add_argument('--thresh_file', required=True, type=str,
            help='pre-built day-of-year percentile threshold file (as used as input to add_mhw_mcs)')
    parser_ts_mhw.add_argument('--out_dir', required=True, type=str,
            help='directory to write one figure per site into')
    parser_ts_mhw.add_argument('--today', required=True, type=str,
            help='the hindcast/forecast transition date, in format "YYYY-MM-DD"')
    parser_ts_mhw.add_argument('--level', required=False, type=parse_int, default=-1,
            help='sigma level to plot: -1 (default) is the surface, 0 is the bottom')
    parser_ts_mhw.add_argument('--depth_name', required=False, type=str, default='Surface',
            help="label for the level, used in the titles and file names e.g. 'Surface' or 'Bottom'")
    def plot_timeseries_mhw_handler(args):
        plot_timeseries_mhw(args.ts_file, args.clim_file, args.thresh_file,
                            args.out_dir, args.today,
                            level=args.level, depth_name=args.depth_name)
    parser_ts_mhw.set_defaults(func=plot_timeseries_mhw_handler)

    # ------------------------------
    # plot_timeseries_stratification
    # ------------------------------
    parser_ts_strat = subparsers.add_parser('plot_timeseries_stratification',
            help='Plot per-site stratification time-series (bottom minus target-depth density)')
    parser_ts_strat.add_argument('--ts_file', required=True, type=str,
            help='time-series file written by get_ts_multivar, containing the stratification variable')
    parser_ts_strat.add_argument('--out_dir', required=True, type=str,
            help='directory to write one figure per site into')
    parser_ts_strat.add_argument('--today', required=True, type=str,
            help='the hindcast/forecast transition date, in format "YYYY-MM-DD"')
    parser_ts_strat.add_argument('--target_depth', required=False, type=float, default=None,
            help='the depth (in metres, positive down) the stratification was computed against, as passed to add_stratification. Only used to label the y axis')
    def plot_timeseries_stratification_handler(args):
        plot_timeseries_stratification(args.ts_file, args.out_dir, args.today,
                                       target_depth=args.target_depth)
    parser_ts_strat.set_defaults(func=plot_timeseries_stratification_handler)

    # ----------------
    # make_clm_inter
    # ----------------
    parser_make_clm_inter = subparsers.add_parser('make_clm_inter',
            help='Make monthly ocean climatology files for CROCO interannual runs')
    parser_make_clm_inter.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the monthly OGCM files')
    parser_make_clm_inter.add_argument('--output_dir', required=True, type=str,
            help='Path to where the forcing files will be saved. This directory also needs a crocotools_param.py file')
    parser_make_clm_inter.add_argument('--month_start', required=True, type=str, 
            help='first month in the interannual run in format "YYYY-MM"')
    parser_make_clm_inter.add_argument('--month_end', required=True, type=str,
            help='last month in the interannual run in format "YYYY-MM"')
    parser_make_clm_inter.add_argument('--Yorig', required=True, type=int,
            help='the Yorig value used in setting up the CROCO model')
    def make_clm_inter_handler(args):
        sys.path.append(args.output_dir)
        import crocotools_param as params

        month_now = datetime.strptime(args.month_start+'-01','%Y-%m-%d')
        month_end = datetime.strptime(args.month_end+'-01','%Y-%m-%d')

        while month_now <= month_end:

            print('working on '+month_now.strftime('%Y-%m'))

            # define a list of 3 input files - last month, this month and next month
            # This is to ensure that we always have clm data for the start and end of the CROCO run for this month
            month_prev = month_now - timedelta(days=10) # an arbitrary date in the previous month
            month_next = month_now + timedelta(days=32) # 32 days ensures we get to the next month
            fname_month_prev = os.path.join(args.input_dir, month_prev.strftime(params.input_file_fmt))
            fname_month_now = os.path.join(args.input_dir, month_now.strftime(params.input_file_fmt))
            fname_month_next = os.path.join(args.input_dir, month_next.strftime(params.input_file_fmt))
            if not os.path.exists(fname_month_prev):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_prev)
            if not os.path.exists(fname_month_now):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_now)
            if not os.path.exists(fname_month_next):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_next)
            input_file=[fname_month_prev,
                        fname_month_now,
                        fname_month_next]

            # define ini_date and end_date using a 1 day buffer either side of the month
            ini_date = month_now - timedelta(days=1)
            day_end = calendar.monthrange(month_now.year,month_now.month)[1]
            end_date = datetime(month_now.year,month_now.month,day_end) + timedelta(days=2)

            # make the clm file for this month
            fname_out = params.clim_prefix + month_now.strftime('_Y%YM%m.nc')
            make_clm(input_file,args.output_dir,ini_date,end_date,args.Yorig,fname_out)

            month_now=datetime(month_next.year, month_next.month, 1) # set month_now to the first day of the next month
        
    parser_make_clm_inter.set_defaults(func=make_clm_inter_handler)
 
    # ----------------
    # croco_srf_2_ww3
    # ----------------
    parser_croco_srf_2_ww3 = subparsers.add_parser('croco_srf_2_ww3',
            help='Convert CROCO surface output to WW3-compatible current and water level netCDF files for use with ww3_prnc in ASIS mode')
    parser_croco_srf_2_ww3.add_argument('--fname', required=True, type=str,
            help='input CROCO surface filename - can include wildcards (*) to process multiple files in a single command')
    parser_croco_srf_2_ww3.add_argument('--grdname', required=False, type=str, default=None,
            help='optional CROCO grid file (if grid vars are not in your output files)')
    parser_croco_srf_2_ww3.add_argument('--dir_out', required=True, type=str,
            help='output directory for WW3-compatible current and water level files')
    parser_croco_srf_2_ww3.add_argument('--Yorig', type=parse_int, default=None,
            help='Origin year for CROCO time. Output time units will be "days since Yorig-01-01 00:00:00"')
    def croco_srf_2_ww3_handler(args):
        from glob import glob as globfn
        fname_in = args.fname
        if type(fname_in) == str:
            if fname_in.find('*') > 0:
                fname_in = sorted(globfn(fname_in))
            else:
                fname_in = [fname_in]
        for f in fname_in:
            croco_srf_2_ww3(f, grdname=args.grdname, dir_out=args.dir_out, Yorig=args.Yorig)
    parser_croco_srf_2_ww3.set_defaults(func=croco_srf_2_ww3_handler)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
