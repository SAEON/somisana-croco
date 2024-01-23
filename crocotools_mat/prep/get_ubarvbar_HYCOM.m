function [ubar,vbar]=...
         get_ubarvbar_HYCOM(hisfile,gridfile,tindex,vlevel,skp,npts,scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% calculate ubar and vbar fileds from 3D HYCOM u and v data
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
if vlevel==0
u=get_hslice(hisfile,gridfile,'ubar',tindex,vlevel,'u');
v=get_hslice(hisfile,gridfile,'vbar',tindex,vlevel,'v');
else
u=get_hslice(hisfile,gridfile,'u',tindex,vlevel,'u');
v=get_hslice(hisfile,gridfile,'v',tindex,vlevel,'v');
end
[u,v,lon_spd,lat_spd,mask_spd]=uv_vec2rho(u,v,lon,lat,angle,mask,skp,npts);
u=u*scale;
v=v*scale;
return
