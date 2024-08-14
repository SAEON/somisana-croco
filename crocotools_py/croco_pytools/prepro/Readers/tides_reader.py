import numpy as np
import sys

def lookvar(input):
    '''
    Dictionary for tides files. How to fill it:
    First, indicate the different grid names for elev and current (u,v) fields
    Then if you data are Amp_Phase type follow the style:
      H_a   Elevation amplitude
      H_p   Elevation phase
      U_a   Eastward velocity amplitude
      U_p   Eastward velocity phase
      V_a   Northward velocity amplitude
      V_p   Northward velocity phase
    If your data are Re_Im type use:
      H_r   Elevation real part
      H_i   Elevation imaginary part
      U_r   Eastward velocity real part
      U_i   Eastward velocity imaginary part
      V_r   Northward velocity real part
      V_i   Northward velocity imaginary part
    '''

    if input == 'tpxo9':
        dico={ 'lonr':'lon_z','lonu':'lon_u','lonv':'lon_v',\
               'latr':'lat_z','latu':'lat_u','latv':'lat_v',\
               'H_r':'hRe',\
               'H_i':'hIm',\
               'U_r':'uRe',\
               'U_i':'uIm',\
               'V_r':'vRe',\
               'V_i':'vIm',\
               'topor':'h','topou':'hu','topov':'hv'\
              }

    elif input == 'tpxo9_lowres':
        dico={ 'lonr':'lon_z','lonu':'lon_u','lonv':'lon_v',\
               'latr':'lat_z','latu':'lat_u','latv':'lat_v',\
               'H_r':'hRe',\
               'H_i':'hIm',\
               'U_r':'URe',\
               'U_i':'UIm',\
               'V_r':'VRe',\
               'V_i':'VIm',\
               'topor':'hz','topou':'hu','topov':'hv'\
              }


    elif input == 'tpxo7_croco':
        dico={ 'lonr':'lon_r','lonu':'lon_u','lonv':'lon_r',\
               'latr':'lat_r','latu':'lat_r','latv':'lat_v',\
               'H_r':'ssh_r',\
               'H_i':'ssh_i',\
               'U_r':'u_r',\
               'U_i':'u_i',\
               'V_r':'v_r',\
               'V_i':'v_i',\
               'topor':'h','topou':'h','topov':'h'\
              }

    elif input == 'fes2014':
        dico={ 'lonr':'lon','lonu':'lon','lonv':'lon',\
               'latr':'lat','latu':'lat','latv':'lat',\
               'H_a':'amplitude',\
               'H_p':'phase',\
               'U_a':'Ua',\
               'U_p':'Ug',\
               'V_a':'Va',\
               'V_p':'Vg'\
              }

    else:
        sys.exit(print('No \'%s\' dico available in Modules/tides_readers.py. Please add it an run the script' % input))

    return dico
