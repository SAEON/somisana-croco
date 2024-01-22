%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a CROCO grid file
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  This is as per make_grid.m in croco_tools, only call to add_topo
%  changed to add_topo_scatter, where bathy is gridded from scatter data
%  e.g. electronic charts
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
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
%  Contributions of P. Marchesiello (IRD) and X. Capet (UCLA)
%
%  Updated    Aug-2006 by Pierrick Penven
%  Updated    24-Oct-2006 by Pierrick Penven (mask correction)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
crocotools_param
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
warning off
isoctave=exist('octave_config_info');
%
% Title
%
if level>0
    grdname=[grdname,'.',num2str(level)];
end
disp(' ')
disp([' Making the grid: ',grdname])
disp(' ')
disp([' Title: ',CROCO_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])
%
% Get the Longitude
%
lonr=(lonmin:dl:lonmax);
%
% Get the latitude for an isotropic grid
%
i=1;
latr(i)=latmin;
while latr(i)<=latmax
  i=i+1;
  latr(i)=latr(i-1)+dl*cos(latr(i-1)*pi/180);
end
[Lonr,Latr]=meshgrid(lonr,latr);
[Lonu,Lonv,Lonp]=rho2uvp(Lonr); 
[Latu,Latv,Latp]=rho2uvp(Latr);
%
% Create the grid file
%
disp(' ')
disp(' Create the grid file...')
[M,L]=size(Latp);
disp([' LLm = ',num2str(L-1)])
disp([' MMm = ',num2str(M-1)])
create_grid(L,M,grdname,CROCO_title)
%
% Fill the grid file
%
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'lat_u'}(:)=Latu;
nc{'lon_u'}(:)=Lonu;
nc{'lat_v'}(:)=Latv;
nc{'lon_v'}(:)=Lonv;
nc{'lat_rho'}(:)=Latr;
nc{'lon_rho'}(:)=Lonr;
nc{'lat_psi'}(:)=Latp;
nc{'lon_psi'}(:)=Lonp;
close(nc)
%
%  Compute the metrics
%
disp(' ')
disp(' Compute the metrics...')
[pm,pn,dndx,dmde]=get_metrics(grdname);
xr=0.*pm;
yr=xr;
for i=1:L
  xr(:,i+1)=xr(:,i)+2./(pm(:,i+1)+pm(:,i));
end
for j=1:M
  yr(j+1,:)=yr(j,:)+2./(pn(j+1,:)+pn(j,:));
end
[xu,xv,xp]=rho2uvp(xr);
[yu,yv,yp]=rho2uvp(yr);
dx=1./pm;
dy=1./pn;
dxmax=max(max(dx/1000));
dxmin=min(min(dx/1000));
dymax=max(max(dy/1000));
dymin=min(min(dy/1000));
disp(' ')
disp([' Min dx=',num2str(dxmin),' km - Max dx=',num2str(dxmax),' km'])
disp([' Min dy=',num2str(dymin),' km - Max dy=',num2str(dymax),' km'])
%
%  Angle between XI-axis and the direction
%  to the EAST at RHO-points [radians].
%
angle=get_angle(Latu,Lonu);
%
%  Coriolis parameter
%
f=4*pi*sin(pi*Latr/180)*366.25/(24*3600*365.25);
%
% Fill the grid file
%
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'pm'}(:)=pm;
nc{'pn'}(:)=pn;
nc{'dndx'}(:)=dndx;
nc{'dmde'}(:)=dmde;
nc{'x_u'}(:)=xu;
nc{'y_u'}(:)=yu;
nc{'x_v'}(:)=xv;
nc{'y_v'}(:)=yv;
nc{'x_rho'}(:)=xr;
nc{'y_rho'}(:)=yr;
nc{'x_psi'}(:)=xp;
nc{'y_psi'}(:)=yp;
nc{'angle'}(:)=angle;
nc{'f'}(:)=f;
nc{'spherical'}(:)='T';
close(nc);
%
%
%  Add topography from topofile
%
disp(' ')
disp(' Add topography...')
h=add_topo_scatter(grdname,topofile,hmin,hmax);
%
% Compute the mask
%
if exist('grd_mask_rho.mat', 'file')
    maskr_str=load('grd_mask_rho.mat');
    maskr=maskr_str.mask_rho;
else 
    maskr=h>0;
end
maskr=process_mask(maskr);
[masku,maskv,maskp]=uvp_mask(maskr);
%
%  Write it down
%
nc=netcdf(grdname,'write');
nc{'h'}(:)=h;
nc{'mask_u'}(:)=masku;
nc{'mask_v'}(:)=maskv;
nc{'mask_psi'}(:)=maskp;
nc{'mask_rho'}(:)=maskr;
close(nc);

if (isoctave == 0)
%
% Create the coastline
%
if ~isempty(coastfileplot)
  make_coast(grdname,coastfileplot);
end
%
r=input('Do you want to use editmask ? y,[n]','s');
if strcmp(r,'y')
  disp(' Editmask:')
  disp(' Edit manually the land mask.')
  disp(' Press enter when finished.')
  disp(' ')
  if ~isempty(coastfileplot)
    editmask(grdname,coastfilemask)
  else
    editmask(grdname)
  end
  r=input(' Finished with edit mask ? [press enter when finished]','s');
end
%
close all
end % isoctave
%
%  Smooth the topography
%
nc=netcdf(grdname,'write');
h=nc{'h'}(:);
mask_rho=nc{'mask_rho'}(:);
%
h=smoothgrid(h,mask_rho,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);
%
%  Write it down
%
disp(' ')
disp(' Write it down...')
nc{'h'}(:)=h;
close(nc);
%
% save the mask
save mask_rho.mat mask_rho

% make a plot
%
% if (isoctave == 0)
% if makeplot==1
%   disp(' ')
%   disp(' Do a plot...')
%   themask=ones(size(maskr));
%   themask(maskr==0)=NaN; 
%   domaxis=[min(min(Lonr)) max(max(Lonr)) min(min(Latr)) max(max(Latr))];
%   colaxis=[min(min(h)) max(max(h))];
%   fixcolorbar([0.25 0.05 0.5 0.03],colaxis,...
%               'Topography',10)
%   width=1;
%   height=0.8;
%   subplot('position',[0. 0.14 width height])
%   m_proj('mercator',...
%          'lon',[domaxis(1) domaxis(2)],...
%          'lat',[domaxis(3) domaxis(4)]);
%   m_pcolor(Lonr,Latr,h.*themask);
%   shading flat
%   caxis(colaxis)
%   hold on
%   [C1,h1]=m_contour(Lonr,Latr,h,[hmin 100 200 500 1000 2000 4000],'k');
%   clabel(C1,h1,'LabelSpacing',1000,'Rotation',0,'Color','r')
%   if ~isempty(coastfileplot)
%     m_usercoast(coastfileplot,'color','r');
%     %m_usercoast(coastfileplot,'speckle','color','r');
%   else
%     m_gshhs_l('color','r');
%     m_gshhs_l('speckle','color','r');
%   end
%   m_grid('box','fancy',...
%          'xtick',5,'ytick',5,'tickdir','out',...
%          'fontsize',7);
%   hold off
% end
% warning on
% end
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

