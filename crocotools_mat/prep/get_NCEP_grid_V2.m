function [imin,imax,jmin,jmax,lon,lat]=...
    get_NCEP_grid_V2(fname,lonmin,lonmax,latmin,latmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts

dl=3;

lonmin=lonmin-dl;
lonmax=lonmax+dl;
latmin=latmin-dl;
latmax=latmax+dl;
%
% Get the global horizontal grid
%
nc=netcdf(fname);
lon=nc{'lon'}(:);
lat=nc{'lat'}(:);
close(nc)
%
% Get a subgrid
%
j=find(lat>=latmin & lat<=latmax);
%
i=find(lon>=lonmin & lon<=lonmax);
%
lon=lon(i);
lat=lat(j);
jmin=min(j);
jmax=max(j);
imin=min(i);
imax=max(i);
%
