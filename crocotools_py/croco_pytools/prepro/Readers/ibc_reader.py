import numpy as np
import sys

def lookvar(input):

    if input == 'mercator':
        dico={ 'depth':'depth',\
               'lonr':'longitude','lonu':'longitude','lonv':'longitude',\
               'latr':'latitude','latu':'latitude','latv':'latitude',\
               'ssh':'zos',\
               'temp':'thetao',\
               'salt':'so',\
               'u': 'uo',\
               'v': 'vo',\
               'time': 'time',\
               'time_dim':'time'\
             }

    elif input == 'eccov4':
        dico={ 'depth':'Z',\
               'lonr':'longitude','lonu':'longitude','lonv':'longitude',\
               'latr':'latitude','latu':'latitude','latv':'latitude',\
               'ssh':'ETAN',\
               'temp': 'THETA',\
               'salt': 'SALT',\
               'u': 'EVEL',\
               'v': 'NVEL',\
               'time':'time',\
               'time_dim':'time'\
             }

    elif input == 'soda':
        dico={ 'depth':'st_ocean',\
               'lonr':'xt_ocean','lonu':'xu_ocean','lonv':'xu_ocean',\
               'latr':'yt_ocean','latu':'yu_ocean','latv':'yu_ocean',\
               'ssh':'ssh',\
               'temp': 'temp',\
               'salt': 'salt',\
               'u': 'u',\
               'v': 'v',\
               'time':'time',\
               'time_dim':'time'\
             }
    elif input == 'mercator_croco':
        dico={ 'depth':'depth',\
               'lonr':'lonT','lonu':'lonU','lonv':'lonV',\
               'latr':'latT','latu':'latU','latv':'latV',\
               'ssh':'ssh',\
               'temp':'temp',\
               'salt':'salt',\
               'u': 'u',\
               'v': 'v',\
               'time': 'time',\
               'time_dim':'time'\
             }
    else:
        sys.exit(print('No \'%s\' dico available in Modules/inputs_readers.py. Please add it an run the script' % input))

    return dico



