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
from crocotools_py.preprocess import make_tides,reformat_gfs_atm,reformat_saws_atm,make_ini,make_bry
from crocotools_py.postprocess import (get_ts_multivar, compute_anomaly,
                                       get_ds, handle_time,
                                       create_mhw_output_netcdf, process_single_level)
from crocotools_py.plotting import plot as crocplot
from crocotools_py.regridding import regrid_tier1, regrid_tier2, regrid_tier3 
 
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
                      cmap = 'Spectral_r',
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
    def regrid_tier2_handler(args):
        regrid_tier2(args.fname, args.dir_out, args.grdname, args.Yorig, args.doi_link, depths = args.depths)
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
    def regrid_tier3_handler(args):
        regrid_tier3(args.fname, args.dir_out, args.Yorig, args.doi_link, spacing=args.spacing)
    parser_regrid_tier3.set_defaults(func=regrid_tier3_handler)
    
    # ----------------
    # get_ts_multivar
    # ----------------
    parser_get_ts_multivar = subparsers.add_parser('get_ts_multivar', help='extract a time-series or a profile through time from a croco file')
    parser_get_ts_multivar.add_argument('--fname', required=True, type=str, help='input CROCO filename')
    parser_get_ts_multivar.add_argument('--lon', required=True, type=float, help='Longitude of data extraction')
    parser_get_ts_multivar.add_argument('--lat', required=True, type=float, help='Latitude of data extraction')
    parser_get_ts_multivar.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be in seconds since Yorig-01-01')
    parser_get_ts_multivar.add_argument('--vars', type=parse_list, 
                        default=['temp', 'salt'],
                        help='optional list of CROCO variable names')
    parser_get_ts_multivar.add_argument('--depths',required=False, type=parse_list,
                        default=[0,-5,-10,-20,-50,-100,-200,-500,-1000,-99999],
                        help='Depths for time-series extraction (see get_ts_multivar() for description of input)')
    parser_get_ts_multivar.add_argument('--fname_out', required=True, help='output filename')
    def get_ts_multivar_handler(args):
        get_ts_multivar(args.fname, args.lon, args.lat, Yorig=args.Yorig, 
               vars = args.vars, 
               depths = args.depths,
               write_nc=True, # default behaviour in the cli is to write a file
               fname_nc=args.fname_out)
    parser_get_ts_multivar.set_defaults(func=get_ts_multivar_handler)
    
    # ----------------
    # compute_anomaly
    # ----------------
    parser_compute_anomaly=subparsers.add_parser('compute_anomaly', help='Compute anomalies by subtracting monthly climatology from high-frequency CROCO output.')
    parser_compute_anomaly.add_argument('--fname_clim', required=True, type=str, help='input climatology file')
    parser_compute_anomaly.add_argument('--fname_in', required=True, type=str, help='input high frequency (forecast) file')
    parser_compute_anomaly.add_argument('--fname_out', required=True, type=str, help='output directory')
    parser_compute_anomaly.add_argument('--Yorig', type=parse_int, 
                        default=2000, 
                        help='Origin year used in setting up CROCO time i.e. CROCO time will be in seconds since Yorig-01-01')
    parser_compute_anomaly.add_argument('--varlist', type=parse_list_str, 
                        default=['temp','u','v', 'salt','zeta'],
                        help='optional list of CROCO variable names')
    parser_compute_anomaly.add_argument('--use_constant_clim', type=parse_bool, 
                       default='False',
                       help='If True, subtract a constant value from fname_in, computed as the climatology at the midpoint of the time axis of fname_in')
    def compute_anomaly_handler(args):
        compute_anomaly(args.fname_clim, args.fname_in, args.fname_out, 
                        args.Yorig, varlist=args.varlist,
                        use_constant_clim=args.use_constant_clim)
    parser_compute_anomaly.set_defaults(func=compute_anomaly_handler)
 
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
    # detect_mhw_forecast
    # ----------------------
    parser_detect_mhw = subparsers.add_parser(
        'detect_mhw_forecast',
        help=(
            'Detect Marine Heatwave (MHW) and Marine Cold Spell (MCS) events '
            'in a CROCO forecast file using a pre-built 4D climatology. '
            'Writes a single NetCDF with signed categories: '
            '+1..+4 = MHW, 0 = no event, -1..-4 = MCS.'))
    parser_detect_mhw.add_argument(
        '--temp_file', required=True, type=str,
        help='Path to the CROCO forecast temperature file (or glob pattern).')
    parser_detect_mhw.add_argument(
        '--clim_file', required=True, type=str,
        help=(
            'Path to the pre-built 4D climatology NetCDF '
            '(Climatology_4D_Unified.nc). Must contain variables '
            'climatology, threshold_90 and threshold_10 with dimensions '
            '(day_of_year, s_rho, eta_rho, xi_rho).'))
    parser_detect_mhw.add_argument(
        '--fname_out', required=True, type=str,
        help='Full path and filename for the output categories NetCDF file.')
    parser_detect_mhw.add_argument(
        '--temp_var', required=False, type=str, default='temp',
        help='Name of the temperature variable in temp_file (default: temp).')
    parser_detect_mhw.add_argument(
        '--Yorig', required=False, type=parse_int, default=2000,
        help=(
            'Reference year for the CROCO time axis in the forecast file '
            'i.e. time is in seconds since Yorig-01-01 (default: 2000).'))
    parser_detect_mhw.add_argument(
        '--batch_size', required=False, type=parse_int, default=5,
        help=(
            'Number of eta_rho rows processed at once. '
            'Reduce if you run out of memory (default: 5).'))

    parser_plot_mhw = subparsers.add_parser(
        'plot_mhw_forecast',
        help='Generate operational MetOcean plots (timeseries, flags, and GIFs) from forecast outputs.')
    parser_plot_mhw.add_argument('--forecast_file', required=True, type=str, help='Path to the raw CROCO forecast file.')
    parser_plot_mhw.add_argument('--cat_file', required=True, type=str, help='Path to the processed Categories NetCDF.')
    parser_plot_mhw.add_argument('--clim_file', required=True, type=str, help='Path to the Climatology NetCDF.')
    parser_plot_mhw.add_argument('--out_dir', required=True, type=str, help='Output directory for the plots.')
    parser_plot_mhw.add_argument('--start_date', required=True, type=str, help='Forecast start date (YYYY-MM-DD)')
    parser_plot_mhw.add_argument('--end_date', required=True, type=str, help='Forecast end date (YYYY-MM-DD)')
    parser_plot_mhw.add_argument('--Yorig', required=False, type=parse_int, default=2000, help='Reference year (default: 2000)')
    
    def plot_mhw_forecast_handler(args):
        # We import here so we don't load huge plotting libraries if running a different cli command
        from crocotools_py.plotting import plot_operational_mhw_mcs
        plot_operational_mhw_mcs(
            args.forecast_file, args.cat_file, args.clim_file, 
            args.out_dir, args.start_date, args.end_date, args.Yorig
        )
    parser_plot_mhw.set_defaults(func=plot_mhw_forecast_handler)
    def detect_mhw_forecast_handler(args):
        import os
        import gc
        import numpy as np
        import pandas as pd
        import xarray as xr
        from pathlib import Path
 
        os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
 
        out_file = Path(args.fname_out)
        out_file.parent.mkdir(parents=True, exist_ok=True)
 
        # -- climatology (lazy, one level loaded per call inside process_single_level)
        print(f'Opening climatology: {args.clim_file}')
        ds_clim    = xr.open_dataset(args.clim_file)
        num_levels = ds_clim.sizes['s_rho']
        n_eta      = ds_clim.sizes['eta_rho']
        n_xi       = ds_clim.sizes['xi_rho']
        doy_values = ds_clim['day_of_year'].values
 
        # -- temperature: get_ds + handle_time avoids the depth computation
        #    that get_var() triggers, which would allocate ~128 GB
        print(f'Loading temperature: {args.temp_file}')
        ds_temp = get_ds(args.temp_file, args.temp_var)
        ds_temp = handle_time(ds_temp, Yorig=args.Yorig)
 
        # Build daily time axis from raw time bounds.
        # Resampling happens per-batch inside process_single_level.
        raw_times    = ds_temp.time.values
        first_day    = pd.Timestamp(raw_times[0]).normalize()
        last_day     = pd.Timestamp(raw_times[-1]).normalize()
        target_dates = pd.date_range(start=first_day, end=last_day, freq='1D')
        T       = len(target_dates)
        t_dates = np.array([d.toordinal() for d in target_dates], dtype=int)
        ds_dummy_daily = xr.Dataset(coords={'time': target_dates})
 
        # -- initialise output
        if out_file.exists():
            out_file.unlink()
 
        nc_out = create_mhw_output_netcdf(
            output_file   = str(out_file),
            n_time_daily  = T,
            n_levels      = num_levels,
            n_eta         = n_eta,
            n_xi          = n_xi,
            ds_temp_daily = ds_dummy_daily,
            ds_temp       = ds_temp,
            ds_clim       = ds_clim,
            mode_name     = 'MHW_MCS',
        )
        cat_var = nc_out.variables['category']
 
        # -- process each level: MHW first, MCS fills no-event days
        try:
            for k in range(num_levels - 1, -1, -1):
                process_single_level(k, num_levels, ds_temp, ds_clim,
                                     args.temp_var, doy_values,
                                     False, t_dates, args.batch_size, nc_out)
                mhw_cats = cat_var[:, k, :, :][:]
 
                cat_var[:, k, :, :] = np.zeros_like(mhw_cats)
                process_single_level(k, num_levels, ds_temp, ds_clim,
                                     args.temp_var, doy_values,
                                     True, t_dates, args.batch_size, nc_out)
                mcs_cats = cat_var[:, k, :, :][:]
 
                combined = mhw_cats.copy()
                combined[combined == 0] = mcs_cats[combined == 0]
                cat_var[:, k, :, :] = combined
 
                nc_out.sync()
                gc.collect()
        finally:
            nc_out.sync()
            nc_out.close()
 
        ds_clim.close()
        size_mb = out_file.stat().st_size / (1024 ** 2)
        print(f'Done: {out_file}  ({size_mb:.1f} MB)')
 
    parser_detect_mhw.set_defaults(func=detect_mhw_forecast_handler)
 
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")
 
if __name__ == "__main__":
    main()
 