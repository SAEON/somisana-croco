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
from crocotools_py.postprocess import get_ts_multivar
from crocotools_py.plotting import plot as crocplot
from crocotools_py.regridding import regrid_tier1, regrid_tier2, regrid_tier3 
from download.cmems import download_glorys, download_mercator
from download.gfs import download_gfs_atm
from download.hycom import download_hycom
from crocotools_py.regridding_cfc import regrid1_cf_compliant,regrid2_cf_compliant,regrid3_cf_compliant

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

def parse_bool(s: str) -> bool:
    try:
        return {'true':True, 'false':False}[s.lower()]
    except KeyError:
        raise argparse.ArgumentTypeError(f'expect true/false, got: {s}')

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-croco repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...

    # ----------------
    # download_glorys
    # ----------------
    parser_download_glorys = subparsers.add_parser('download_glorys', 
            help='Download month by month of daily MERCATOR 1/12 deg reanalysis data (GLORYS) from CMEMS')
    parser_download_glorys.add_argument('--usrname', required=True, type=str, help='Copernicus username')
    parser_download_glorys.add_argument('--passwd', required=True, help='Copernicus password')
    parser_download_glorys.add_argument('--domain', type=parse_list, 
                        default=[23, 34, -37, -31],
                        help='comma separated list of domain extent to download i.e. "lon0,lon1,lat0,lat1"')
    parser_download_glorys.add_argument('--start_date', required=True, type=parse_datetime, 
                        help='start time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_glorys.add_argument('--end_date', required=True, type=parse_datetime, 
                        help='end time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_glorys.add_argument('--varList', type=parse_list, 
                        default=['so', 'thetao', 'zos', 'uo', 'vo'],
                        help='comma separated list of variables to download e.g. "so,thetao,zos,uo,vo"')
    parser_download_glorys.add_argument('--depths', type=parse_list, 
                        default=[0.493, 5727.918],
                        help='comma separated list of depth extent to download (positive down). For all depths use "0.493,5727.918"')
    parser_download_glorys.add_argument('--outputDir', required=True, help='Directory to save files')
    def download_glorys_handler(args):
        download_glorys(args.usrname, args.passwd, args.domain, args.start_date,args.end_date,args.varList, args.depths, args.outputDir)
    parser_download_glorys.set_defaults(func=download_glorys_handler)
    
    # -------------------
    # download_mercator
    # -------------------
    parser_download_mercator = subparsers.add_parser('download_mercator', 
            help='Download a subset of daily MERCATOR 1/12 deg analysis data from CMEMS')
    parser_download_mercator.add_argument('--usrname', required=True, type=str, help='Copernicus username')
    parser_download_mercator.add_argument('--passwd', required=True, help='Copernicus password')
    parser_download_mercator.add_argument('--domain', type=parse_list, 
                        default=[10, 25, -40, -25],
                        help='comma separated list of domain extent to download i.e. "lon0,lon1,lat0,lat1"')
    parser_download_mercator.add_argument('--run_date', required=True, type=parse_datetime, 
                        help='start time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_mercator.add_argument('--hdays', required=True, type=float,
                        default=5.,
                        help='hindcast days i.e before run_date')
    parser_download_mercator.add_argument('--fdays', required=True, type=float,
                        default=5.,
                        help='forecast days i.e before run_date')
    parser_download_mercator.add_argument('--outputDir', required=True, help='Directory to save files') 
    def download_mercator_handler(args):
        download_mercator(args.usrname, args.passwd, args.domain, args.run_date,args.hdays, args.fdays,args.outputDir)
    parser_download_mercator.set_defaults(func=download_mercator_handler)
    
    # ------------------
    # download_gfs_atm
    # ------------------
    parser_download_gfs_atm = subparsers.add_parser('download_gfs_atm', 
            help='Download a subset of hourly 0.25 deg gfs atmospheric data from NOAA')
    parser_download_gfs_atm.add_argument('--domain', type=parse_list, 
                        default=[10, 25, -40, -25],
                        help='comma separated list of domain extent to download i.e. "lon0,lon1,lat0,lat1"')
    parser_download_gfs_atm.add_argument('--run_date', required=True, type=parse_datetime, 
                        help='start time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_gfs_atm.add_argument('--hdays', required=True, type=float,
                        default=5.,
                        help='hindcast days i.e before run_date')
    parser_download_gfs_atm.add_argument('--fdays', required=True, type=float,
                        default=5.,
                        help='forecast days i.e before run_date')
    parser_download_gfs_atm.add_argument('--outputDir', required=True, help='Directory to save files')
    def download_gfs_atm_handler(args):
        download_gfs_atm(args.domain, args.run_date, args.hdays, args.fdays, args.outputDir)
    parser_download_gfs_atm.set_defaults(func=download_gfs_atm_handler)
    
    # -------------------
    # download_hycom
    # -------------------
    parser_download_hycom = subparsers.add_parser('download_hycom', 
            help='Download a subset of  HYCOM analysis data using xarray OpenDAP')
    parser_download_hycom.add_argument('--outDir', required=True, help='Directory to save files')
    parser_download_hycom.add_argument('--domain', required = False, type=parse_list,
            default=[10.0, 25.0, -40.0, -20.0],
            help='comma separated list of domain extent to download i.e. "lon0,lon1,lat0,lat1"')
    parser_download_hycom.add_argument('--depths', required = False, type=parse_list,
            default=[0,2,4,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,200,250,300,350,400,500,600,700,800,900,1000,1250,1500,2000,2500,3000,4000,5000],
            help='comma separated list of depths to download i.e. "[0,2,4,...,3000,4000,5000]"')
    parser_download_hycom.add_argument('--variables', required = False, type=parse_list,
            default=['surf_el','salinity','water_temp','water_v','water_u'],
            help='List of strings of variable names found in HYCOM to download i.e.["surf_el","salinity","water_temp","water_u","water_v"].')
    parser_download_hycom.add_argument('--run_date', required=False, type=parse_datetime,
            default=None,
            help='start time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_hycom.add_argument('--hdays', required=False, type=float, 
            default=5.,
            help='hindcast days i.e before run_date')
    parser_download_hycom.add_argument('--fdays', required=False, type=float, 
            default=5.,
            help='forecast days i.e before run_date')
    parser_download_hycom.add_argument('--cleanDir', required=False, type=parse_bool, 
            default=True,
            help='Clean the directory after merging the files')
    parser_download_hycom.add_argument('--parallel', required=False, type=parse_bool, 
            default=True,
            help='Type of download. If parallel, then the download occurs in parellel. If parallel is false, then the download occurs in series. ')
    def download_hycom_handler(args):
        download_hycom(args.outDir,args.domain, args.depths, args.variables, args.run_date, args.hdays, args.fdays,  args.cleanDir, args.parallel)
    parser_download_hycom.set_defaults(func=download_hycom_handler)
    
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
    parser_crocplot.add_argument('--gif_out', required=True, type=str, help='the output gif filename')
    parser_crocplot.add_argument('--level', required=False, default=None, type=parse_int, help='level to plot. If >=0, then a sigma level is plotted. If <0 then a z level (in m) is plotted. Default behaviour will plot the surface layer')
    parser_crocplot.add_argument('--ticks', required=False, type=parse_list,
                         default=[12,13,14,15,16,17,18,19,20,21,22], 
                         help='contour ticks to use in plotting the variable')
    parser_crocplot.add_argument('--cbar_label', required=False, default='temperature ($\degree$C)', type=str, help='the label used for the colorbar')
    parser_crocplot.add_argument('--isobaths', required=False, type=parse_list,
                         default=[100,500],
                         help='the isobaths to add to the figure')
    parser_crocplot.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    def crocplot_handler(args):
        # a lot of the crocplot inputs are hard coded below but could be made configurable in future
        # there are also other potential optional inputs to this function which we aren't specifying here
        # so this is a work in progress...
        crocplot(args.fname, var=args.var,
                      grdname=None, # could make this configurable in the cli args
                      tstep=0,
                      tstep_end=None,
                      level=args.level,
                      ticks = args.ticks,
                      cmap = 'Spectral_r',
                      extents = None,
                      ref_date = args.ref_date,
                      cbar_label=args.cbar_label,
                      add_vectors = True,
                      skip_time = 1,
                      isobaths=args.isobaths,
                      gif_out=args.gif_out,
                      write_gif = True
                      )
    parser_crocplot.set_defaults(func=crocplot_handler)
    
    # --------------
    # regrid_tier1
    # --------------
    parser_regrid_tier1 = subparsers.add_parser('regrid_tier1', 
            help='tier 1 regridding of a raw CROCO output file: regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid -> rotates u/v to be east/north components instead of grid-aligned components -> adds a depth variable providing the depths of each sigma level at each time-step')
    parser_regrid_tier1.add_argument('--fname', required=True, type=str, help='input native CROCO filename')
    parser_regrid_tier1.add_argument('--fname_out', required=True, help='tier 1 output filename')
    parser_regrid_tier1.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    def regrid_tier1_handler(args):
        regrid_tier1(args.fname, args.fname_out, args.ref_date)
    parser_regrid_tier1.set_defaults(func=regrid_tier1_handler)
    
    # --------------
    # regrid_tier2
    # --------------
    parser_regrid_tier2 = subparsers.add_parser('regrid_tier2', 
            help='tier 2 regridding of a CROCO output: takes the output of regrid-tier1 as input and regrids the sigma levels to constant z levels, including the surface and bottom layers -> output variables are the same as tier 1, only depths is now a dimension with the user specified values')
    parser_regrid_tier2.add_argument('--fname', required=True, type=str, help='input regridded tier1 filename')
    parser_regrid_tier2.add_argument('--fname_out', required=True, help='tier 2 output filename')
    parser_regrid_tier2.add_argument('--depths', required=False, type=parse_list,
                         default=[0,-5,-10,-20,-50,-100,-200,-500,-1000,-99999],  
                         help='list of depths to extract (in metres, negative down). A value of 0 denotes the surface and a value of -99999 denotes the bottom layer)')
    def regrid_tier2_handler(args):
        regrid_tier2(args.fname, args.fname_out, depths = args.depths)
    parser_regrid_tier2.set_defaults(func=regrid_tier2_handler)
    
    # --------------
    # regrid_tier3
    # --------------
    parser_regrid_tier3 = subparsers.add_parser('regrid_tier3', 
            help='tier 3 regridding of a CROCO output: takes the output of regrid-tier3 as input and regrids the horizontal grid to be regular with a specified grid spacing. Output variables are the same as tier 1 and 2, only horizontal grid is now rectilinear with hz dimensions of longitude,latitude i.e. horizontal grid is no longer curvilinear. The extents of the rectilinear grid are automatically determined using the curvilinear grid extents.')    
    parser_regrid_tier3.add_argument('--fname', required=True, type=str, help='input regridded tier2 filename')
    parser_regrid_tier3.add_argument('--fname_out', required=True, help='tier 3 output filename')
    parser_regrid_tier3.add_argument('--spacing', type=str,
                         default='0.01', 
                         help='constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output')
    def regrid_tier3_handler(args):
        regrid_tier3(args.fname, args.fname_out, spacing = args.spacing)
    parser_regrid_tier3.set_defaults(func=regrid_tier3_handler)
    
    # ----------------
    # get_ts_multivar
    # ----------------
    parser_get_ts_multivar = subparsers.add_parser('get_ts_multivar', help='extract a time-series or a profile through time from a croco file')
    parser_get_ts_multivar.add_argument('--fname', required=True, type=str, help='input CROCO filename')
    parser_get_ts_multivar.add_argument('--lon', required=True, type=float, help='Longitude of data extraction')
    parser_get_ts_multivar.add_argument('--lat', required=True, type=float, help='Latitude of data extraction')
    parser_get_ts_multivar.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    parser_get_ts_multivar.add_argument('--vars', type=parse_list, 
                        default=['temp', 'salt'],
                        help='optional list of CROCO variable names')
    parser_get_ts_multivar.add_argument('--depths', type=parse_list,
                        default=[0,-5,-10,-20,-50,-100,-200,-500,-1000,-99999],
                        help='Depths for time-series extraction (see get_ts_multivar() for description of input)')
    parser_get_ts_multivar.add_argument('--fname_out', required=True, help='output filename')
    def get_ts_multivar_handler(args):
        get_ts_multivar(args.fname, args.lon, args.lat, args.ref_date, 
               vars = args.vars, 
               depths = args.depths,
               write_nc=True, # default behaviour in the cli is to write a file
               fname_nc=args.fname_out)
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
            
            fname_out = params.croco_prefix + month_now.strftime('_Y%YM%m.nc')
            
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
        
        fname_out = params.ini_prefix + ini_date.strftime('_Y%YM%m.nc')
        
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
    
    # ----------------
    # Regrid Tier 1 CF-Compliant
    # ----------------
    parser_regrid1_cfc = subparsers.add_parser('regrid1_cfc',
            help='Tier 1 regridding of a raw CROCO model output file/s: regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid -> rotates u/v to be east/north components instead of grid-aligned components -> adds a depth variable providing the depths of each sigma level at each time-step -> saves outputs in a CF-Compliant netCDF that can be published')
    parser_regrid1_cfc.add_argument('--fname', required=True, type=str,nargs='+', 
            help='Native CROCO filename. Can be filename (fname = /path/to/file.nc), list of files (fname = [file1,file2,file2]) or wildcard (fname = /path/to/*.nc).')
    parser_regrid1_cfc.add_argument('--info_dir', required=True, 
            help='Path to directory containing the configuration information (crocotools_param.py and croco_cf_compliance.py).')
    parser_regrid1_cfc.add_argument('--out_dir', required=False, default=None, 
            help='Save directory. If out_dir is None, then a new directory is made where fname is.')
    def regrid1_cfc_handler(args):
        regrid1_cf_compliant(args.fname, args.info_dir, args.out_dir)
    parser_regrid1_cfc.set_defaults(func=regrid1_cfc_handler)
    
    # ----------------
    # Regrid Tier 2 CF-Compliant
    # ----------------
    parser_regrid2_cfc = subparsers.add_parser('regrid2_cfc',
            help='Tier 2 regridding of a CROCO output: takes the output of regrid-tier1 as input and regrids the sigma levels to constant z-levels, including the surface and bottom layers -> output variables are the same as tier 1, only depths is now a dimension with the user specified values. Outputs are stored in CF-Compliant netCDF files that can be pubplished.')
    parser_regrid2_cfc.add_argument('--fname', required=True, type=str,nargs='+',
            help='Native CROCO filename. Can be filename (fname = /path/to/file.nc), list of files (fname = [file1,file2,file2]) or wildcard (fname = /path/to/*.nc).')
    parser_regrid2_cfc.add_argument('--info_dir', required=True,
            help='Path to directory containing the configuration information (crocotools_param.py and croco_cf_compliance.py).')
    parser_regrid2_cfc.add_argument('--depths',required=False,default=None,type=parse_list,
            help='list of depths to extract (in metres, negative down). A value of 0 denotes the surface and a value of -99999 denotes the bottom layer). Usually this is set at: depths = [0,-5,-10,-20,-50,-100,-200,-500,-1000,-99999]')
    parser_regrid2_cfc.add_argument('--out_dir', required=False,default=None,
            help='Save directory. If out_dir is None, then a new directory is made where fname is.')
    def regrid2_cfc_handler(args):
        regrid2_cf_compliant(args.fname, args.info_dir, args.depths, args.out_dir)
    parser_regrid2_cfc.set_defaults(func=regrid2_cfc_handler)
    
    # ----------------
    # Regrid Tier 3 CF-Compliant
    # ----------------
    parser_regrid3_cfc = subparsers.add_parser('regrid3_cfc',
            help='Tier 3 regridding of a CROCO output: takes the output of regrid-tier2 as input and regrids the curvilinear grid to rectilinear -> output variables are the same as tier 1 and 2, only now longitude and lititude are 1-Dimensional arrays -> Outputs are stored in CF-Compliant netCDF files that can be pubplished -> Outputs are mostly used for visualisation and web display and not for science due to the interpolation.')
    parser_regrid3_cfc.add_argument('--fname', required=True, type=str,nargs='+',
            help='Input filename from tier 2 regridding. Can be filename (fname = /path/to/file.nc), list of files (fname = [file1,file2,file2]) or wildcard (fname = /path/to/*.nc).')
    parser_regrid3_cfc.add_argument('--info_dir', required=True,
            help='Path to directory that contains the configuration information (crocotools_param.py and croco_cf_compliance.py).')
    parser_regrid3_cfc.add_argument('--spacing',  required=False, type=float, default=None,
            help='constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output')
    parser_regrid3_cfc.add_argument('--out_dir', required=False,default=None,
            help='Save directory. If out_dir is None, then a new directory is made where fname is.')
    def regrid3_cfc_handler(args):
        regrid3_cf_compliant(args.fname, args.info_dir, args.spacing, args.out_dir)
    parser_regrid3_cfc.set_defaults(func=regrid3_cfc_handler)
    


    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
