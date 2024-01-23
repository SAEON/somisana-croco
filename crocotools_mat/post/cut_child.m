function var_out=...
         cut_child(var,lon,lat,lon_1,lat_1,npts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% set points inside domain to NaN where child domain data exists
%
% var = 2D matrix with variable data
% lon = parent lon
% lat = parent lat
% lon_1 = child lon
% lat_1 = child lat
% npts = number of grid points around border of child to exclude from the
%        cut. Format is [S E N W]. Useful to include some overlap to avoid 
%        white space between parent and child domains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6
    npts=[1 1 1 1];
end
% define the corners of the child domain 
% this code is from bounddomain.m
[Mp,Lp]=size(lon_1);
% cut grid cells in to avoid white space around child when plotting
% parent and child together
imin=1+npts(4);
imax=Lp-npts(2);
jmin=1+npts(1);
jmax=Mp-npts(3);
xsquare=cat(1,lon_1(jmin:jmax,imin),  ...
                lon_1(jmax,imin:imax)' ,...
                lon_1(jmax:-1:jmin,imax),...
                lon_1(jmin,imax:-1:imin)' );
ysquare=cat(1,lat_1(jmin:jmax,imin),  ...
                lat_1(jmax,imin:imax)' ,...
                lat_1(jmax:-1:jmin,imax),...
                lat_1(jmin,imax:-1:imin)' );

%determine which parent domain points are inside the child domain
in = inpolygon(lon,lat,xsquare,ysquare); % contains 1's where points are inside
out = ones(size(in));
out = out - in; % contains 1's where points are outside

var_out = var.*out;
var_out(var_out==0)=NaN; % sets all points inside child domain to missing

return
