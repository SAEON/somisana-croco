%--------------------------------------------------------------
%
%  Modify child grid topography such that it matches the interpolated
%  parent topography at the boundaries.
%  
%  This script is for use with make_roms2roms.
%
%   (c) 2007 Jeroen Molemaker, UCLA
%
%  original script modified by GGF to be run with crocotools_param.m
%  this version uses the HYCOM 1/12 degree global model bathymetry
%  as the parent and blends it with the bathy contained in the parent croco
%  grid
%--------------------------------------------------------------
%

clear all
close all
crocotools_param

chdgrd=grdname; % child grid in this case is by definition the level 0 
                % roms grid, as the parent is an external model
pargrd=global_topofile; % parent grid 
                % in this case is the hycom grid

% End user-defined----------------------------------------------
%

% Get minimal parent subgrid bounds
lonc = ncread(chdgrd, 'lon_rho')';
latc = ncread(chdgrd, 'lat_rho')';

dl = 0.1;
lon0 = min(min(lonc)-dl);
lon1 = max(max(lonc))+dl;
lat0 = min(min(latc))-dl;
lat1 = max(max(latc))+dl;

% plot(lonp(jmin:jmax,imin:imax),latp(jmin:jmax,imin:imax),'.k')
% hold on;plot(lonc,latc,'.r');hold off
% drawnow

% Get topography data from childgrid
hc   = ncread(chdgrd, 'h')';
mask = ncread(chdgrd, 'mask_rho')';
[Mc,Lc]=size(hc);

% Get parent grid and squeeze minimal subgrid
nc=netcdf(pargrd,'r');
x=nc{'longitude'}(:);
y=nc{'latitude'}(:);
%  get a subgrid
j=find(y>=lat0 & y<=lat1);
i=find(x>=lon0 & x<=lon1);
x=x(i);
y=y(j);
[lonp,latp,fp,pmp,pnp]=get_tpx_grid(x,y);
hp=nc{'deptho'}(j,i);
hp_missing = ncreadatt(pargrd,'deptho','_FillValue');
hp(hp==hp_missing)=NaN;
close(nc)
maskp = hp>0;

% Get interpolation coefficient to go to (lonc,latc).
[elem,coef] = get_tri_coef(lonp,latp,lonc,latc,maskp);
%% parent grid topo at child locations
hpi = sum(coef.*hp(elem),3);

% very NB step we need to include here (according to GGF anyway) is that we
% need to replace values in hpi which are NaN in parent interpolation, but
% not in the child- just use child values in this case
nan_indx=isnan(hpi);
hpi(nan_indx)=hc(nan_indx);
% now adjust to limit at hmax
hpi(hpi>hmax)=hmax;

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
