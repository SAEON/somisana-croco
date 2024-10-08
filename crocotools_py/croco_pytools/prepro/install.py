#!/usr/bin python3
# -*- coding: utf-8 -*-
"""
Main croco_pytools install execuable script
"""

import os
import subprocess
from __init__ import (
    ENV_MOD,
    ENV_CONDA,
)

install_venv_conda = input(
        "Do you want to install conda environment? [y,[n]]: "
                          )

fortran_compilation = input( 
        "Do you want to compile fortran tools? [y,[n]]:"
                           )

if install_venv_conda.lower() in ['y','yes']:
    venv_conda = ''.join((ENV_CONDA,'/environment_tools.yml'))
    print(f"Install conda env from {venv_conda}")
    os.chdir(ENV_CONDA)

    mamb_install = input(
        "Do you want use mamba to create the environnemnt? [y,[n]]: "
                         )
    if mamb_install.lower() in ['y','yes']:
        try:
            subprocess.run(f"mamba env create -q -f {venv_conda}",
                       shell=True,
                       check=True)
        except:
            try:
                print('Installation with mamba failled. Trying with conda')
                subprocess.run(f"conda env create -q -f {venv_conda}",
                           shell=True,
                           check=True)
            except:
                print(f"Could not find {venv_conda}")
    else:
        try:
            subprocess.run(f"conda env create -q -f {venv_conda}",
                       shell=True,
                       check=True)
        except:
            print(f"Could not find {venv_conda}")
 
else:
    print('No installation done')



if fortran_compilation.lower() in ['y','yes']:
    print('Compiling fortran tools')
    env_fine = ['croco_pyenv','envcroco']
    for i,ev in enumerate(env_fine):
        try:
            subprocess.run(f"conda list --name {ev}",
                            shell=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)
            print(f"Found conda environment '{ev}' to compile fortran tools")
            env_to_use = ev
            break
        except:
            env_to_use = ''
            continue

    if env_to_use in env_fine:
        os.chdir(ENV_MOD+'/tools_fort_routines')
#        subprocess.run(f"conda run -n {env_to_use} python compilation_fortran_tools.py",
        subprocess.run(["make", "clean"], check=True)
        log = subprocess.run(f"conda run -n {env_to_use} make",
                       shell=True,
                       check=True)
        print(log)
