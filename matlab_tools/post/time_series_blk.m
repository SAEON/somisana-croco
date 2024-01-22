function var=time_series_blk(blkfile,gridfile,lon0,lat0,vname,coef)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  Get a timeseries from a roms bulk file (uesful for checking input)
%
%  bashed on time_series.m as part of roms_tools
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Defaults values
%
if nargin < 1
  error('You must specify a file name')
end
if nargin < 2
  gridfile=blkfile;
  disp(['Default grid name: ',gridfile])
end
if nargin < 3
  lon0=[];
end
if nargin < 4
  lat0=[];
end
if nargin < 5
  vname='wspd';
  disp(['Default variable to plot: ',vname])
end
if nargin < 6
  coef=1;
  disp(['Default coef: ',num2str(coef)])
end
%
% Get default values
%
if isempty(gridfile)
  gridfile=blkfile;
end
if strcmp(vname,'uwnd') || strcmp(vname,'sustr')
  [lat,lon,mask]=read_latlonmask(gridfile,'u');
elseif strcmp(vname,'vwnd') || strcmp(vname,'svstr')
  [lat,lon,mask]=read_latlonmask(gridfile,'v');
else
  [lat,lon,mask]=read_latlonmask(gridfile,'r');
end
if isempty(lon0) || isempty(lat0)
  lat0=mean(mean(lat));
  lon0=mean(mean(lon));
end
%
% Find j,i indices for the timeseries
%
disp(['lon0 = ',num2str(lon0),' - lat0 = ',num2str(lat0)])
[J,I]=find((lat(1:end-1,1:end-1)<=lat0 & lat(2:end,2:end)>lat0 &...
            lon(2:end,1:end-1)<=lon0 & lon(1:end-1,2:end)>lon0)==1);
if isempty(I) ||  isempty(J)
  disp('Warning: timeseries location not found')
  [M,L]=size(lon);
  I=round(L/2);
  J=round(M/2);
end
% hack edit in case two points are found
I=I(1);
J=J(1);
disp(['I = ',int2str(I),' J = ',int2str(J)])
lon1=lon(J,I);
lat1=lat(J,I);
disp(['lon1 = ',num2str(lon1),' - lat1 = ',num2str(lat1)])
%
% get the vertical levels
%
nc=netcdf(blkfile);
var=coef*squeeze(nc{vname}(:,J,I));
close(nc)

return

