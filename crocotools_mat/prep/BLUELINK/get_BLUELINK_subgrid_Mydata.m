function [i1min,i1max,i2min,i2max,i3min,i3max,jrange,krange,lonT,latT,lonUV,latUV,depth]=...
    get_BLUELINK_subgrid_Mydata(bluelinkfile_t,bluelinkfile_u,lonmin,lonmax,latmin,latmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  Get the indices for a BLUELINK subgrid 
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
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dl=1;
lonmin=lonmin-dl;
lonmax=lonmax+dl;
latmin=latmin-dl;
latmax=latmax+dl;
%
% Get the grid
%
nc_t=netcdf(bluelinkfile_t);
latT=nc_t{'yt_ocean'}(:);
lonT=nc_t{'xt_ocean'}(:);
depth=nc_t{'st_ocean'}(:);
nc_t(close);
nc_u=netcdf(bluelinkfile_u);
latUV=nc_u{'yu_ocean'}(:);
lonUV=nc_u{'xu_ocean'}(:);
nc_u(close);
%
% Get a subgrid
%
% 1 Longitude: take care of greenwitch
%
i1=find(lonT-360>=lonmin & lonT-360<=lonmax);
i2=find(lonT>=lonmin & lonT<=lonmax);
i3=find(lonT+360>=lonmin & lonT+360<=lonmax);
%
lonT=cat(1,lonT(i1)-360,lonT(i2),lonT(i3)+360);
% GGF comment - we'll just use the rho grid to get the subset indices - 
% these indices will be fine for subsetting the u and v grids too (no need
% to over-complicate with separate indices for the separate grids as we're 
% only using these indices to subset roughly, so long as the subset covers
% our croco grid)
lonUV=cat(1,lonUV(i1)-360,lonUV(i2),lonUV(i3)+360);
%
if ~isempty(i1)
  i1min=min(i1);
  i1max=max(i1);
else
  i1min=[];
  i1max=[];
end
if ~isempty(i2)
  i2min=min(i2);
  i2max=max(i2);
else
  i2min=[];
  i2max=[];
end
if ~isempty(i3)
  i3min=min(i3);
  i3max=max(i3);
else
  i3min=[];
  i3max=[];
end
%
% 2 Latitude
%
j=find(latT>=latmin & latT<=latmax);
latT=latT(j);
latUV=latUV(j);
jmin=min(j);
jmax=max(j);
jrange=jmin:jmax;
%
% 3 Depth
%
k=length(depth);
krange=1:k;
%
return
