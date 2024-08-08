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
from datetime import datetime
from crocotools_py.preprocess import reformat_gfs_atm
from crocotools_py.postprocess import get_ts_multivar
from crocotools_py.regridding import regrid_tier1, regrid_tier2, regrid_tier3 
from download.cmems import download_glorys, download_mercator
from download.gfs import download_gfs_atm
from download.hycom import download_hycom

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYY-MM-DD HH:MM:SS'.")

def parse_list(value):
    return [x.strip() for x in value.split(',')]

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
                        default=[23, 34, -37, -31],
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
                        default=[23, 34, -37, -31],
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
    parser_download_hycom.add_argument('--domain', required = False, type=parse_list,
            default=[23.0, 34.0, -37.0, -31.0],
            help='comma separated list of domain extent to download i.e. "lon0,lon1,lat0,lat1"')
    parser_download_hycom.add_argument('--depths', required = False, type=parse_list,
            default=[0,2,4,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100,125,150,200,250,300,350,400,500,600,700,800,900,1000,1250,1500,2000,2500,3000,4000,5000],
            help='comma separated list of depths to download i.e. "[0,2,4,...,3000,4000,5000]"')
    parser_download_hycom.add_argument('--variables', required = False, type=parse_list,
            default=['surf_el','salinity','water_temp','water_v','water_u'],
            help='List of strings of variable names found in HYCOM to download i.e.["surf_el","salinity","water_temp","water_u","water_v"].')
    parser_download_hycom.add_argument('--run_date', required=True, type=parse_datetime,
            help='start time in format "YYYY-MM-DD HH:MM:SS"')
    parser_download_hycom.add_argument('--hdays', required=False, type=float, default=5.,
            help='hindcast days i.e before run_date')
    parser_download_hycom.add_argument('--fdays', required=False, type=float, default=5.,
            help='forecast days i.e before run_date')
    parser_download_hycom.add_argument('--outDir', required=True, help='Directory to save files')
    parser_download_hycom.add_argument('--cleanDir', required=False, type=bool, default=True,
            help='Clean the directory after merging the files')
    parser_download_hycom.add_argument('--parallel', type=bool, required=False, default=True,
            help='Type of download. If parallel, then the download occurs in parellel. If parallel is false, then the download occurs in series. ')
    def download_hycom_handler(args):
        download_hycom(args.domain, args.depths, args.variables, args.run_date, args.hdays, args.fdays, args.outDir, args.cleanDir, args.parallel)
    parser_download_hycom.set_defaults(func=download_hycom_handler)
    
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
    parser_regrid_tier2.add_argument('--depths', type=parse_list,
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

    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    
