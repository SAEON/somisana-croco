This repo is part of the [SOMISANA](https://somisana.ac.za/) initiative, and used for our [CROCO](https://www.croco-ocean.org/) related model development.

Directories in the repository: 
- `crocotools_mat`:    some matlab code which is used in conjunction with the official croco-tools to support some of SOMISANA's model development, particularly the pre-processing. 
- `crocotools_py`:     python functions for some postprocessing, validation and plotting CROCO model output
- `download`:          python functions for downloading data used for forcing our CROCO configurations
- `configs`:           configurations used for SOMISANA's hindcast and forecast simulations (see README's in the sub-directories for further details)
- `.github/workflows`: github workflows for running our forecast models operationally on a server set up for this purpose on MIMS (`run_ops.yml` is the `main` workflow)

This repo is largely a redesign of [this one](https://github.com/SAEON/somisana) (big shout out to Zach Smith and Matt Carr for their work on this). The new repo is more 'model-centric', and we don't include anything website related, which was integrated into the model development in the previous design.

Please refer to the [wiki](https://github.com/SAEON/somisana-croco/wiki) for more information.

