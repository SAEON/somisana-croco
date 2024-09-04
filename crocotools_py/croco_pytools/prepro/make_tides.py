__author__ = 'Mathieu Le Corre'
__email__  = 'mathieu.le.corre@shom.fr'
__date__   = '2022-11'
__license__='GPL3'

'''
===========================================================================
Further Information:  
  http://www.croco-ocean.org
  
This file is part of CROCOTOOLS

Create a CROCO forcing file for tides
In the current state the script can handle:
    - TPXO (raw data)
    - FES2014

To add a new dataset you just have to go in Readers/tides_readers.py and
create a dico. Follow information given inside the file to fill it properly 

The periodicity of a dataset (if it is -180/180 or 0/360) is handle by the
script. Just download the dataset and it should be fine

The script works as follow:
    - reads croco_grd
    - reads input data and restricts the area of spatial coverage
    - creates croco_frc.nc
    - Loop on choosen waves
    - Writes data in netcdf

===========================================================================
'''

#--- Dependencies ---------------------------------------------------------

import netCDF4 as netcdf
import cftime
import numpy as np
import glob
import sys
from os import path
sys.path.append("./Modules/")
sys.path.append("./Readers/")
import interp_tools
import croco_class as Croco
import tides_class as Inp

#--- USER CHANGES ---------------------------------------------------------

# Dates
# Origin year
Yorig = 2000 # 1900 if TIDES_MAS defined in cppdef.h
# Initial date
Yini, Mini, Dini = 2013, 1, 1


# Input data information and formating
# Note: if you are using a tpxo dataset please be sure to have somewhere in 
#       inputdata 'tpxo'. This will allow the code to use the OTIS (TPXO is obtained with it)
#       convention a-b*i for complex.
#       Also, if you have already preprocess TPXO, make sure that you have correct units and 
#       u,v are in m/s and not m**2/s
inputdata = 'tpxo7_croco' # Input data dictionnary as defined in the Readers/tides_reader.py
input_dir = '../../DATASETS_CROCOTOOLS/TPXO7/'
input_file = 'TPXO7.nc' # Leave empty if you have multiple files
input_type = 'Re_Im' # Format of the input data 'Amp_phase'or 'Re_Im'
multi_files  = False # Set to True if several input files
if multi_files: 
    waves_separated = True # Set to True if input files waves are separated
    elev_file = 'h_<tides>_tpxo9_atlas_30_v5.nc' # elevation file names. if wave_separated put <tides> where wave name is found
    u_file = 'u_<tides>_tpxo9_atlas_30_v5.nc' # eastward currents file names. if wave_separated put <tides> where wave name is found
    v_file = 'u_<tides>_tpxo9_atlas_30_v5.nc' # northward currents file names. if wave_separated put <tides> where wave name is found

# CROCO grid informations
croco_dir = '../../CROCO_FILES/'
croco_grd = 'croco_grd.nc'

# Tide file informations
croco_filename = 'croco_frc.nc'
tides = ['M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm']

cur = True # Set to True if you to compute currents
pot = False # Set to True if you to compute potiential tides

# Nodal correction
Correction_ssh = True
Correction_uv = True

#--- END USER CHANGES -----------------------------------------------------

#--- START MAIN SCRIPT ----------------------------------------------------

#read tides and periods
tides_param=np.loadtxt("./Modules/tides.txt",skiprows=3,comments='^',usecols=(0,4),dtype=str)
tides_names=tides_param[:,0].tolist()
tides_periods=tides_param[:,1].tolist()
tides_names=[x.lower() for x in tides_names]
tides_names=np.array(tides_names)

# Load croco_grd
sigma_params = dict(theta_s=0, theta_b=0, N=1, hc=1) # no needs as 2d
crocogrd = Croco.CROCO_grd(''.join((croco_dir, croco_grd)), sigma_params)

# --- Load input (restricted to croco_grd) ----------------------------

if multi_files:
    if waves_separated:
        input_file_ssh=[]
        input_file_u=[]
        input_file_v=[]
        for inp in tides:
            if path.isfile(input_dir+elev_file.replace('<tides>',inp)):
                input_file_ssh+=[input_dir+elev_file.replace('<tides>',inp)]
            elif path.isfile(input_dir+elev_file.replace('<tides>',inp.lower())):
                input_file_ssh+=[input_dir+elev_file.replace('<tides>',inp.lower())]
            else:
                sys.exit('Elevation file %s for wave %s is missing' % (input_dir+elev_file.replace('<tides>',inp), inp))
           
            if cur:
                if path.isfile(input_dir+u_file.replace('<tides>',inp)):
                    input_file_u+=[input_dir+u_file.replace('<tides>',inp)]
                elif path.isfile(input_dir+u_file.replace('<tides>',inp.lower())):
                    input_file_u+=[input_dir+u_file.replace('<tides>',inp.lower())]
                else:
                    sys.exit('Eastward current file for wave %s is missing' % inp)

                if path.isfile(input_dir+v_file.replace('<tides>',inp)):
                    input_file_v+=[input_dir+v_file.replace('<tides>',inp)]
                elif path.isfile(input_dir+v_file.replace('<tides>',inp.lower())):
                    input_file_v+=[input_dir+v_file.replace('<tides>',inp.lower())]
                else:
                    sys.exit('Northward current file for wave %s is missing' % inp)
        
        input_file_ssh = list(input_file_ssh)
        if cur:
            input_file_u = list(input_file_u)
            input_file_v = list(input_file_v)
        else: 
            input_file_u = None
            input_file_v = None
    else:
        input_file_ssh=list(input_dir+elev_file)
        if cur:
            input_file_u = list(input_dir+u_file)
            input_file_v = list(input_dir+v_file)
        else:
            input_file_u = None
            input_file_v = None
else:
    input_file_ssh=list([input_dir+input_file])
    if cur:
        input_file_u=list([input_dir+input_file])
        input_file_v=list([input_dir+input_file])
    else:
        input_file_u=None
        input_file_v=None

inpdat=Inp.getdata(inputdata,input_file_ssh,crocogrd,input_type,tides,input_file_u,input_file_v)

# --- Create the initial file -----------------------------------------

Croco.CROCO.create_tide_nc(None,''.join((croco_dir+croco_filename)),crocogrd,cur=cur,pot=pot)

if Correction_ssh or Correction_uv:
    date=cftime.datetime(Yini,Mini,Dini)
    date_orig=cftime.datetime(Yorig,1,1)

todo=['H']

if cur:
    todo+=['cur']
if pot:
    todo+=['pot']
    coslat2=np.cos(np.deg2rad(crocogrd.lat))**2
    sin2lat=np.sin(2.*np.deg2rad(crocogrd.lat))

# --- Start loop on waves --------------------------------------------
nc=netcdf.Dataset(croco_dir+croco_filename, 'a')

if Correction_ssh or Correction_uv:
    nc.Nodal_Correction=''.join(('Origin time is ',str(date_orig)))
else:
    nc.Nodal_Correction='No nodal correction'

for i,tide in enumerate(tides) :
    print('\nProcessing *%s* wave' %(tide))
    print('-----------------------')
    index=np.argwhere(tides_names==tide.lower())
    if (len(index)>0):
        print("  tides %s is in the list"%(tide))
        # get period
        index=index[0][0]
        period=float(tides_periods[index])
        print("  Period of the wave %s is %f"%(tide,period))
        
        nc.variables['tide_period'][i]=period
    if multi_files:
        # In this case waves had been concatenated in the order they appear in the list
        tndx=i
    else:
        # read ntime/periods dimension and find the closest wave
        tndx=np.argwhere(abs(inpdat.ntides-period)<1e-4)
        if len(tndx)==0:
            sys.exit('  Did not find wave %s in input file' % tide)
        else:
            tndx=tndx[0]

    [pf,pu,mkB]=inpdat.egbert_correction(tide,date)
    # For phase shift time should be in seconds relatively Jan 1 1992
    # As mkB is the wave phase at this date
    t0 = cftime.date2num(date_orig,'seconds since 1992-01-01:00:00:00')
    if Correction_ssh or Correction_uv:      
        correc_amp   = pf
        correc_phase = mkB+np.deg2rad(t0/(period*10)) +         pu
        #              |--- phase at origin time ---|  |nodal cor at ini time|
    else:
        correc_amp=1
        correc_phase= mkB+np.deg2rad(t0/(period*10))

    # --- Start loop on var ------------------------------------------
    for vars in todo:
        # get data
        if vars == 'H':
            print('\n  Processing tidal elevation')
            print('  -------------------------')
            (tide_complex,NzGood) = interp_tools.interp_tides(inpdat,vars,-1,crocogrd,tndx,tndx,input_type)
            if Correction_ssh:
                tide_amp=np.ma.abs(tide_complex)*correc_amp
                if 'tpxo' in inputdata :
                    tide_phase=np.mod(np.ma.angle(tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                else:
                    tide_phase=np.mod(np.ma.angle(tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
            else:
                tide_amp=np.ma.abs(tide_complex)
                if 'tpxo' in inputdata :
                    tide_phase=np.mod(np.ma.angle(tide_complex)*-180/np.pi,360)
                else:
                    tide_phase=np.mod(np.ma.angle(tide_complex)*180./np.pi,360)  

            nc.variables['tide_Ephase'][i,:]=tide_phase*crocogrd.maskr
            nc.variables['tide_Eamp'][i,:]=tide_amp*crocogrd.maskr

        #########################
        elif vars == 'cur':
            print('\n  Processing tidal currents')
            print('  -------------------------')
            (u_tide_complex,v_tide_complex,NzGood) = interp_tools.interp_tides(inpdat,vars,-1,crocogrd,tndx,tndx,input_type)
            
            if Correction_uv:
                u_tide_amp=np.ma.abs(u_tide_complex)*correc_amp
                v_tide_amp=np.ma.abs(v_tide_complex)*correc_amp
 
                if 'tpxo' in inputdata:
                    u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                    v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*-180/np.pi-correc_phase*180/np.pi,360)
                else:
                    u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
                    v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*180./np.pi-correc_phase*180/np.pi,360)
            else:
                u_tide_amp=np.ma.abs(u_tide_complex)
                v_tide_amp=np.ma.abs(v_tide_complex)

                if 'tpxo' in inputdata:
                    u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*-180/np.pi,360)
                    v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*-180/np.pi,360)
                else:
                    u_tide_phase=np.mod(np.ma.angle(u_tide_complex)*180./np.pi,360)
                    v_tide_phase=np.mod(np.ma.angle(v_tide_complex)*180./np.pi,360)


            major,eccentricity,inclination,phase=inpdat.ap2ep(u_tide_amp,u_tide_phase,v_tide_amp,v_tide_phase)
  
            nc.variables['tide_Cmin'][i,:,:]=major[:,:]*eccentricity[:,:]*crocogrd.maskr
            nc.variables['tide_Cmax'][i,:,:]=major[:,:]*crocogrd.maskr
            nc.variables['tide_Cangle'][i,:,:]=inclination[:,:]*crocogrd.maskr
            nc.variables['tide_Cphase'][i,:,:]=phase[:,:]*crocogrd.maskr
        #########################
        elif vars == 'pot':
            print('\n  Processing equilibrium tidal potential')
            print('  --------------------------------------')
            try:
                coef=eval(''.join(('inpdat.pot_tide.',tide.lower())))
            except:
                try:
                    # some waves start with a number (ex: 2N2) and python do not like it
                    coef=eval(''.join(('inpdat.pot_tide._',tide.lower())))
                except:
                    print('No potential prameter defined for wave %s' % input_wav)
                    coef=[1 ,0]

            if period<13:  # semidiurnal
                Pamp=correc_amp*coef[0]*coef[1]*coslat2
                Ppha=np.mod(-2*crocogrd.lon-correc_phase*180/np.pi,360)
            elif period<26: # diurnal
                Pamp=correc_amp*coef[0]*coef[1]*sin2lat;
                Ppha=np.mod(-crocogrd.lon-correc_phase*180/np.pi,360)
            else: # long-term
                Pamp=correc_amp*coef[0]*coef[1]*(1-1.5*coslat2);
                Ppha=np.mod(-correc_phase*180/np.pi,360.0)
 
            nc.variables['tide_Pamp'][i,:,:]   = Pamp*crocogrd.maskr
            nc.variables['tide_Pphase'][i,:,:] = Ppha*crocogrd.maskr

nc.close()



