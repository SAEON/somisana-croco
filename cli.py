'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow functions to be run from the docker image for this repo
So this is the entry point for the docker image (see Dockerfile_py).
But it's also handy if you want to execute a python function from inside a bash script
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import argparse
from datetime import datetime
from crocotools_py.postprocess import get_ts, get_ts_uv, get_profile, get_profile_uv
from crocotools_py.regridding import regrid_tier1, regrid_tier2, regrid_tier3 
from download.cmems import download_glorys

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYY-MM-DD HH:MM:SS'.")

def parse_list(value):
    return [x.strip() for x in value.split(',')]

# here we just execute the functions using whatever arguments have been povided by the user
def download_glorys_handler(args):
    download_glorys(args.usrname, 
                                args.passwd, 
                                args.domain, 
                                args.start_date,
                                args.end_date,
                                args.varList, 
                                args.depths, 
                                args.outputDir)

def regrid_tier1_handler(args):
    regrid_tier1(args.fname, args.fname_out, args.ref_date)
    
def regrid_tier2_handler(args):
    regrid_tier2(args.fname, args.fname_out, depths = args.depths)
    
def regrid_tier3_handler(args):
    regrid_tier3(args.fname, args.fname_out, spacing = args.spacing)

def get_ts_handler(args):
    get_ts(args.fname, args.var, args.lon, args.lat, args.ref_date, 
           depth = args.depth,
           write_nc=True, # default behaviour in the cli is to write a file
           fname_nc=args.fname_out)

def get_ts_uv_handler(args):
    get_ts_uv(args.fname, args.lon, args.lat, args.ref_date, 
           depth = args.depth,
           write_nc=True, # default behaviour in the cli is to write a file
           fname_nc=args.fname_out)

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-croco repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...

    # Subparser for download_glorys
    parser_download_glorys = subparsers.add_parser('download_glorys', 
            help='Download month by month of daily MERCATOR 1/12 deg reanalysis data (GLORYS)')
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
    parser_download_glorys.set_defaults(func=download_glorys_handler)

    # Subparser for regrid_tier1
    parser_regrid_tier1 = subparsers.add_parser('regrid_tier1', 
            help='tier 1 regridding of a raw CROCO output file: regrids u/v to the density (rho) grid so all parameters are on the same horizontal grid -> rotates u/v to be east/north components instead of grid-aligned components -> adds a depth variable providing the depths of each sigma level at each time-step')
    parser_regrid_tier1.add_argument('--fname', required=True, type=str, help='input native CROCO filename')
    parser_regrid_tier1.add_argument('--fname_out', required=True, help='tier 1 output filename')
    parser_regrid_tier1.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    parser_regrid_tier1.set_defaults(func=regrid_tier1_handler)
    
    # Subparser for regrid_tier2
    parser_regrid_tier2 = subparsers.add_parser('regrid_tier2', 
            help='tier 2 regridding of a CROCO output: takes the output of regrid-tier1 as input and regrids the sigma levels to constant z levels, including the surface and bottom layers -> output variables are the same as tier 1, only depths is now a dimension with the user specified values')
    parser_regrid_tier2.add_argument('--fname', required=True, type=str, help='input regridded tier1 filename')
    parser_regrid_tier2.add_argument('--fname_out', required=True, help='tier 2 output filename')
    parser_regrid_tier2.add_argument('--depths', type=str,
                         default='0,1,2,5,10,15,20,30,40,50,60,70,100,150,200,500,1000,1500,2000,99999', 
                         help='string with comma separated depth levels (in metres, positive down) to interpolate to. You should include a value of 0 to denote the surface and a value of 99999, which indicates the bottom layer)')
    parser_regrid_tier2.set_defaults(func=regrid_tier2_handler)
    
    # Subparser for regrid_tier3
    parser_regrid_tier3 = subparsers.add_parser('regrid_tier3', 
            help='tier 3 regridding of a CROCO output: takes the output of regrid-tier3 as input and regrids the horizontal grid to be regular with a specified grid spacing. Output variables are the same as tier 1 and 2, only horizontal grid is now rectilinear with hz dimensions of longitude,latitude i.e. horizontal grid is no longer curvilinear. The extents of the rectilinear grid are automatically determined using the curvilinear grid extents.')    
    parser_regrid_tier3.add_argument('--fname', required=True, type=str, help='input regridded tier2 filename')
    parser_regrid_tier3.add_argument('--fname_out', required=True, help='tier 3 output filename')
    parser_regrid_tier3.add_argument('--spacing', type=str,
                         default='0.01', 
                         help='constant horizontal grid spacing (in degrees) to be used for the horizontal interpolation of the output')
    parser_regrid_tier3.set_defaults(func=regrid_tier3_handler)
    
    # Subparser for get_ts
    parser_get_ts = subparsers.add_parser('get_ts', help='extract a time-series from a croco file')
    parser_get_ts.add_argument('--fname', required=True, type=str, help='input native filename')
    parser_get_ts.add_argument('--var', required=True, type=str, help='CROCO variable name')
    parser_get_ts.add_argument('--lon', required=True, type=float, help='Longitude of data extraction')
    parser_get_ts.add_argument('--lat', required=True, type=float, help='Latitude of data extraction')
    parser_get_ts.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    parser_get_ts.add_argument('--depth', required=True, type=int,help='Depth for time-series extraction (see get_ts() for description of input)')
    parser_get_ts.add_argument('--fname_out', required=True, help='output time-series filename')
    parser_get_ts.set_defaults(func=get_ts_handler)
    
    # Subparser for get_ts_uv
    parser_get_ts_uv = subparsers.add_parser('get_ts_uv', help='extract a time-series of u,v current components from a croco file')
    parser_get_ts_uv.add_argument('--fname', required=True, type=str, help='input native filename')
    parser_get_ts_uv.add_argument('--lon', required=True, type=float, help='Longitude of data extraction')
    parser_get_ts_uv.add_argument('--lat', required=True, type=float, help='Latitude of data extraction')
    parser_get_ts_uv.add_argument('--ref_date', type=parse_datetime, 
                        default=datetime(2000,1,1,0,0,0), 
                        help='CROCO reference date in format "YYYY-MM-DD HH:MM:SS"')
    parser_get_ts_uv.add_argument('--depth', required=True, type=int,help='Depth for time-series extraction (see get_ts() for description of input)')
    parser_get_ts_uv.add_argument('--fname_out', required=True, help='output time-series filename')
    parser_get_ts_uv.set_defaults(func=get_ts_uv_handler)

    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
    