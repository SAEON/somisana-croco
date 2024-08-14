Building the environnement
--------------------------

Installation with conda and using the "install" script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: 

  Be sure to have conda install on your computer. If not, refer to
  `Conda website <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ to install it.

.. note:: 
  
  the conda environment contains all the elements needed to compile fortran 
  tools:
  ::

      - gcc>=8.3.0
      - gfortran>=8.3.0
      - netcdf-fortran=4.6.0
    

* Execute the installation script
  ::

    python install.py

* Answer the question on the screen:
  ::

    Do you want to install conda environment? [y,[n]]: y
    Do you want to compile fortran tools? y,[n]: y

If envrionment installation is selected, you will be asked an additional 
question:

::

  Do you want use mamba to create the environnemnt? [y,[n]]:

.. note::

  To use this option you need to have 
  `mamba <https://mamba.readthedocs.io/en/latest/>`_ package. you can install it
  with the command:

  ::
 
    conda install -n base --override-channels -c conda-forge mamba

If no error is raised environment installation and compilation were successful.
Now you can load your conda environment:
::

  conda activate croco_pyenv
 
That's it, you are ready to use the preprocessing tools!
Remind to activate this envrironment each time you want to use the 
preprocessing tools.

Installation with another python package manager (e.g. micromamba) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the ``environment_tools.yml`` file provided in ``env``
::

    micromamba create -f env/environment_tools.yml  -c conda-forge

Activate the environment: 
::

    micromamba activate croco_pyenv

Fortran tools compilation (if not using the install.py)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only if you have not used the ``install.py`` script (first section), or if you have answer ``no`` to the question:
``Do you want to compile fortran tools? y,[n]`` when using ``install.py``, you can compile the fortran tools following this procedure.

Activate your python environment:
::

    conda activate croco_pyenv

Or:
::

    micromamba activate croco_pyenv

.. note:: 
  
  the croco_pyenv environment contains all the elements needed to compile fortran 
  tools:
  ::

      - gcc>=8.3.0
      - gfortran>=8.3.0
      - netcdf-fortran=4.6.0

Launch the Fortran routines compilation:
::

    cd Modules/tools_fort_routines
    make clean
    make

If successfull you should now have this file in ``Modules``:
::

    toolsf.cpython-39-x86_64-linux-gnu.so


Common errors
^^^^^^^^^^^^^

You might face some errors while trying to compile fortran tools. 
Here is a list of what have been already encountered by some users and 
the associated solution.

* ifort can raise wn error while compiling. In ``Modules/tools_fort_routines/Makedefs`` try to add 
  ``--fcompiler=intelem`` in ``FFLAGS``.

* ImportError means you have missing librairies. In your terminal do ``nf-config 
  --flibs`` and check that you have ``-lnetcdff -lnetcdf``. 

