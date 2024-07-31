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
%  Updated    28-dec-2020 by Andres Sepulveda (octave/octcdf)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
more off
warning off
pkg load octcdf

disp(['Add the paths of the different toolboxes'])
% the somisana matlab code is in this repo so we're adding using a relative path
addpath(genpath('../../../../crocotools_mat/misc'));
addpath('../../../../crocotools_mat/prep/');
addpath('../../../../crocotools_mat/prep/GLORYS/');
% point to wherever your official croco_tools is
tools_path='/home/gfearon/code/croco_tools-v1.3.1/';
myutilpath=[tools_path,'UTILITIES/'];

%disp(['Add the paths of the different toolboxes'])
%tools_path='../croco_tools/';
%croco_path='../croco/';
%myutilpath=[tools_path,'UTILITIES/'];
%
% Other software directories
%
addpath([myutilpath,'m_map1.4h'])
addpath([myutilpath,'air_sea'])
addpath([myutilpath,'mask'])
%
% ROMSTOOLS directories
%
addpath([tools_path,'Aforc_CFSR'])
addpath([tools_path,'Aforc_ERA5'])
addpath([tools_path,'Aforc_ECMWF'])
addpath([tools_path,'Aforc_NCEP'])
addpath([tools_path,'Aforc_QuikSCAT'])
addpath([tools_path,'Diagnostic_tools'])
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
