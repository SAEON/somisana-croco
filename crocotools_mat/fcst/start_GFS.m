%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add the paths of the different toolboxes
% 
%  this file is here purely for the case of when we need to run reformat_GFS.m
%  this is done outside the context of any particular croco configuration
%  so we need to have a start file here to define the paths
%  the intention is to use this file when we run our matlab docker image. 
%  So paths are relative to the docker image 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Add the paths of the different toolboxes'])
% paths are relative to the docker image used as part of somisana ops
% adjust for local development
addpath(genpath('/somisana-croco/crocotools_mat/misc')); % path inside our docker image
addpath('/somisana-croco/crocotools_mat/fcst/'); % path inside our docker image
% point to wherever your official croco_tools is
tools_path='/croco_tools-v1.3.1/'; % path inside our docker image
myutilpath=[tools_path,'UTILITIES/'];
%
% Other software directories
% you need to have copied the contents of <DATASETS_CROCOTOOLS>/m_map1.4f/private/* 
% (available for download via the croco website) 
% into <myutilpath>/m_map1.4h/private/
% so that the high resolution coastline is available 
addpath([myutilpath,'m_map1.4h']) 
addpath([myutilpath,'air_sea'])
addpath([myutilpath,'mask'])
%
% CROCOTOOLS directories
%
addpath([tools_path,'Aforc_CFSR'])
addpath([tools_path,'Aforc_ERA5'])
addpath([tools_path,'Aforc_ECMWF'])
addpath([tools_path,'Aforc_NCEP'])
addpath([tools_path,'Aforc_QuikSCAT'])
addpath([tools_path,'Diagnostic_tools'])
addpath([tools_path,'Diagnostic_tools/EKE'])
addpath([tools_path,'Forecast_tools'])
addpath([tools_path,'Nesting_tools'])
addpath([tools_path,'Preprocessing_tools'])
addpath([tools_path,'Preprocessing_tools/Bio'])
addpath([tools_path,'Oforc_OGCM'])
addpath([tools_path,'Tides'])
addpath([tools_path,'Tides/T_TIDE'])
addpath([tools_path,'Visualization_tools'])
addpath([tools_path,'Rivers'])
addpath([tools_path,'Town'])
%
%-------------------------------------------------------
%
addpath([myutilpath,'mexcdf/mexnc'])   % 32 and 64 bits version of mexnc 
%
% - If these directories are already in your matlab native path, 
% you can comment these lines
addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf'])
addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncsource'])
addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/nctype'])
addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncutility'])
%
%-------------------------------------------------------
