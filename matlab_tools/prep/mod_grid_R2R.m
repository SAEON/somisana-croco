%--------------------------------------------------------------
%
%  Modify child grid topography such that it matches the interpolated
%  parent topography at the boundaries.
%  
%  This script is for use with make_croco2croco.
%
%   (c) 2007 Jeroen Molemaker, UCLA
%
%  original script modified by GGF to be run with crocotools_param.m
%--------------------------------------------------------------
%

clear all
close all
crocotools_param

if level>0
    grdname=[grdname,'.',num2str(level)];
end

chdgrd=grdname;
pargrd=grdname_p;

% End user-defined----------------------------------------------
%

% Get minimal parent subgrid bounds
lonc = ncread(chdgrd, 'lon_rho')';
latc = ncread(chdgrd, 'lat_rho')';

lonp = ncread(pargrd, 'lon_rho')';
latp = ncread(pargrd, 'lat_rho')';

lon0 = min(min(lonc));
lon1 = max(max(lonc));
lat0 = min(min(latc));
lat1 = max(max(latc));

g = lonp>=lon0&lonp<=lon1 & latp>=lat0&latp<=lat1;
jmin = min(find(any(g')));
jmax = max(find(any(g')));
imin = min(find(any(g)));
imax = max(find(any(g)));
clear g
lj = length(jmin:jmax);
li = length(imin:imax);

% plot(lonp(jmin:jmax,imin:imax),latp(jmin:jmax,imin:imax),'.k')
% hold on;plot(lonc,latc,'.r');hold off
% drawnow

% Get topography data from childgrid
hc   = ncread(chdgrd, 'h')';
mask = ncread(chdgrd, 'mask_rho')';
[Mc,Lc]=size(hc);

% Get parent grid and squeeze minimal subgrid

% GGF comment- not sure why ncread didn't like the original code (now
% commented out, as seems as per documented use for subsetting??? Had to
% edit it to do the subsetting on the next line

%hp    = ncread(pargrd, 'h', [imin jmin], [li lj])';
hp    = ncread(pargrd, 'h')';
hp = hp(jmin:jmax,imin:imax);
%lonp    = ncread(pargrd, 'lon_rho', [imin jmin], [li lj])';
lonp    = ncread(pargrd, 'lon_rho')';
lonp = lonp(jmin:jmax,imin:imax);
%latp    = ncread(pargrd, 'lat_rho', [imin jmin], [li lj])';
latp    = ncread(pargrd, 'lat_rho')';
latp = latp(jmin:jmax,imin:imax);
%maskp    = ncread(pargrd, 'mask_rho', [imin jmin], [li lj])';
maskp    = ncread(pargrd, 'mask_rho')';
maskp = maskp(jmin:jmax,imin:imax);
lonp(lonp<0) = lonp(lonp<0);

% Get interpolation coefficient to go to (lonc,latc).
[elem,coef] = get_tri_coef(lonp,latp,lonc,latc,maskp);
%% parent grid topo at child locations
hpi = sum(coef.*hp(elem),3);

% very NB step we need to include here (according to GGF anyway) is that we
% need to replace values in hpi which are NaN in parent interpolation, but
% not in the child- just use child values in this case
nan_indx=isnan(hpi);
hpi(nan_indx)=hc(nan_indx);

dist = zeros(Mc,Lc,4);
for i = 1:Mc   %% north south
    for j = 1:Lc    %% east west
        dist(i,j,1) =      i/Mc + (1-obc(1))*1e6; % South
        dist(i,j,2) = (Lc-j)/Lc + (1-obc(2))*1e6; % East
        dist(i,j,3) = (Mc-i)/Mc + (1-obc(3))*1e6; % North
        dist(i,j,4) =      j/Lc + (1-obc(4))*1e6; % West
    end
end
dist = min(dist,[],3); % dist is now matrix covering the child grid with values 
% increasing from zero at the boundaries to 0.5 at the centre

% Feel free to play with this function.
%alpha = 0.5*tanh(100*(dist-0.03))+0.5;  % this one ramps from 0 to 1 over about
                                         % 5% of the model domain lenth
alpha = 0.5*tanh( 50*(dist-0.06))+0.5; % this one ramps from 0 to 1 over about
                                       % 10% of the model domain lenth

hcn = alpha.*hc + (1-alpha).*hpi;

ncwrite(chdgrd,'h',hcn');

%% Visualize the modification
%   sc0 = min(min(hcn));
%   sc1 = max(max(hcn));
%   subplot(2,2,1)
%   pcolor(lonc,latc,hpi);caxis([sc0 sc1]);colorbar;shading flat
%   title('Interpolated Parent Topo')
%   subplot(2,2,2)
%   pcolor(lonc,latc,hcn);caxis([sc0 sc1]);colorbar;shading flat
%   title('Boundary Smoothed Child Topo')
%   subplot(2,2,3)
%   pcolor(lonc,latc,hcn-hpi);colorbar;shading flat
%   title('Difference between Parent and child Topo');
%   subplot(2,2,4)
%   pcolor(lonc,latc,alpha);colorbar;shading flat
%   title('Parent/Child transition function');
