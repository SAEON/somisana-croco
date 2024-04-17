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
% Get the path to the mexcdf (it depends on the architecture)
% Comment  all these lines if you don't want to pass in these tests
!uname -m > .mysystem
fid=fopen('.mysystem');
mysystem=fscanf(fid,'%s');

if ( strcmp(mysystem(end-1:end),'86') )
 mysystem2='32';
elseif ( strcmp(mysystem(end-1:end),'64') )
 mysystem2='64';
end

fclose(fid);
matversion=version('-release');
myversion=str2num(matversion(1:2));
!rm -f .mysystem
disp(['Arch : ',mysystem,' - Matlab version : ',matversion])


if ((myversion > 13)    )
  disp(['Use of mexnc and loaddap in ',mysystem2,' bits.'])
  addpath([myutilpath,'mexcdf/mexnc'])   % 32 and 64 bits version of mexnc
%
% - If these directories are already in your matlab native path,
% you can comment these lines
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncsource'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/nctype'])
  addpath([myutilpath,'mexcdf/netcdf_toolbox/netcdf/ncutility'])
%
% Use of built in opendap libraries (no loaddap) - S. Illig 2015
%
  addpath([tools_path,'Opendap_tools_no_loaddap'])
%
%-------------------------------------------------------
elseif (myversion <= 13)
  disp('Use of mex60 and loaddap in 32 bits.')
  addpath([myutilpath,'mex60'])         % Older/32 bits version of mexcdf

% - If these directories are already in your matlab native path,
% you can comment these lines
% - In this case, if problems with subsrefs.m ans subsasign.m,
% it is because there is a conflict with another native subs.m routines in the
% symbolic native toolbox

  addpath([myutilpath,'netcdf_matlab_60'])
  addpath([myutilpath,'netcdf_matlab_60/nctype'])
  addpath([myutilpath,'netcdf_matlab_60/ncutility'])
%
% Use of loaddap  (older versions of matlab)
%
  addpath([tools_path,'Opendap_tools'])

else
  disp(['Arch : ',mysystem,...
       ' you should provide the paths of your own loaddap and mexcdf directories'])
end
