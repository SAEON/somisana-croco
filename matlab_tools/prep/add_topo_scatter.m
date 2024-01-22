function h=add_topo_scatter(grdname,toponame,hmin,hmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% add a topography to a ROMS grid
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% this script is based on add_topo.m, but in this case
% scatter data saved in a ascii file is interpolated onto the roms grid
%
% Last update Pierrick Penven 8/2006.
%
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
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr 
%
%  Updated    Aug-2006 by Pierrick Penven
%  Updated    2006/10/05 by Pierrick Penven (dl depend of model
%                                           resolution at low resolution)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  read roms grid
%
nc=netcdf(grdname,'r');
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
close(nc);
%
% Get ROMS averaged resolution
%
dx=mean(mean(1./pm));
dy=mean(mean(1./pn));
dx_roms=mean([dx dy]);
disp(['   ROMS resolution : ',num2str(dx_roms/1000,3),' km'])
%
%  open the topo file
%
%  get lon/lat/depth arrays from scatter data
[~,~,fileExt] = fileparts(toponame);
if strcmpi(fileExt, '.csv')
    % Read data using csvread if the file has a .csv extension
    scatterdata = csvread(toponame);
else
    % Read data as space-separated values if the file has a different extension
    scatterdata = dlmread(toponame, ' ');
end

% 
%  interpolate the scatter data onto the roms grid
%
h=griddata(scatterdata(:,1),scatterdata(:,2),scatterdata(:,3),lon,lat,'natural');
% the above interpolation can result in NaN's over the area of interest
% as griddata returns NaN for query points outside of the convex hull (not
% sure what that means). So we just fill the gaps with nearest neighbour
% interpolation which will always give you full coverage
h_nearest=griddata(scatterdata(:,1),scatterdata(:,2),scatterdata(:,3),lon,lat,'nearest');
nan_indx=isnan(h);
h(nan_indx)=h_nearest(nan_indx);
h(h<hmin)=hmin;
h(h>hmax)=hmax;
%
return
