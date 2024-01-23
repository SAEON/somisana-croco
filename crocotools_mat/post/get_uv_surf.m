function [u,v,lon_spd,lat_spd,mask_spd]=...
         get_uv_surf(hisfile,gridfile,tindex,skp,npts,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% put a uv current field in the carthesian frame
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[lat,lon,mask]=read_latlonmask(gridfile,'r');
nc=netcdf(gridfile);
angle=nc{'angle'}(:);
if isempty(angle)
disp('Warning: no angle found in history file')
angle=0*lat;
end
close(nc);
%vlevel=length(ncreadatt(hisfile,'/','sc_r')); % so that we can just use get_hslice below 
%u=get_hslice(hisfile,gridfile,'u',tindex,vlevel,'u');
%v=get_hslice(hisfile,gridfile,'v',tindex,vlevel,'v');
nc=netcdf(hisfile);
u=squeeze(nc{'surf_u'}(tindex,:,:));
v=squeeze(nc{'surf_v'}(tindex,:,:));
close(nc)
[u,v,lon_spd,lat_spd,mask_spd]=uv_vec2rho(u,v,lon,lat,angle,mask,skp,npts);
u=u*scale;
v=v*scale;
return
