function [i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax,lon,lat]=...
    get_CFSR_hourly_subgrid(ifile,lonmin,lonmax,latmin,latmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Get the indices for a CFSR hourly download subgrid 
%
% from get_ECMWF_subgrid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
dl=3;
lonmin=lonmin-dl;
lonmax=lonmax+dl;
latmin=latmin-dl;
latmax=latmax+dl;
%
% Get the global grid
%
nc=netcdf(ifile);
lon=squeeze(nc{'lon'}(:));
lat=squeeze(nc{'lat'}(:));
close(nc);
%
% Get a subgrid
%
%
% 1 Longitude: take care of greenwitch
%
i1=find(lon-360>=lonmin & lon-360<=lonmax);
i2=find(lon>=lonmin & lon<=lonmax);
i3=find(lon+360>=lonmin & lon+360<=lonmax);
%
lon=cat(1,lon(i1)-360,lon(i2),lon(i3)+360);
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
j=find(lat>=latmin & lat<=latmax);
lat=lat(j);
jmin=min(j);
jmax=max(j);
%
% North-South inversion (! after jmin and jmax !)
%
lat=flipdim(lat,1);
%
return
