%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add the paths of the different toolboxes
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Patrick Marchesiello and Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    10-Sep-2006 by Pierrick Penven
%  Updated    22-Sep-2006 by Pierrick Penven (64 bits test)
%  Updated    24-Oct-2006 by Pierrick Penven (mask added)
%  Updated    16-jan-2007 by Pierrick Penven (quikscat added)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Add the paths of the different toolboxes'])
% the somisana matlab code is in this repo so we're adding using a relative path
addpath(genpath('../../../../crocotools_mat/misc'));
addpath('../../../../crocotools_mat/prep/GRID');
% point to wherever your official croco_tools is
tools_path='/home/gfearon/code/croco_tools-v1.3.1/';
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
