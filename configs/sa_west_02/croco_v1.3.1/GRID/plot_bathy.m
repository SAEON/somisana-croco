clear all;
close all;

start

lon_sub=ncread('croco_grd.nc','lon_rho');
lat_sub=ncread('croco_grd.nc','lat_rho');
mask_sub=ncread('croco_grd.nc','mask_rho');
dl=0.5;
lon_submin=min(min(lon_sub))-dl;
lon_submax=max(max(lon_sub))+dl;
lat_submin=min(min(lat_sub))-dl;
lat_submax=max(max(lat_sub))+dl;
h=ncread('croco_grd.nc','h');
% assign masked points to missing
themask=ones(size(mask_sub));
themask(mask_sub==0)=NaN;
h=themask.*h;
h=-h;

% contours
contours=[-4000:500:-500 -300 -200 -150 -100 -50 -40 -30 -20 -10];
c_ticks = -4000:20:0; %This needs to be updated for the new depth range
no_c=length(c_ticks)-1;
c_lims=[min(c_ticks) max(c_ticks)];
%
h(h<min(contours))=min(contours);

Fig = figure();
%resize figure
scale_width=1.3;
scale_height=1;
s=Fig.Position;
s(3)=s(3)*scale_width;
s(4)=s(4)*scale_height;
set(Fig,'Position',s);
%figureFullScreen(Fig);
set(Fig,'color','w');
ha = tight_subplot(1,1,[.05 .05],[.11 .05],[.04 .02]);
axes(ha(1))

% set up figure extents
lonmin=lon_submin;
lonmax=lon_submax;
latmin=lat_submin;
latmax=lat_submax;

domaxis=[lonmin lonmax latmin latmax];

m_proj('mercator','lon',[domaxis(1) domaxis(2)], 'lat',[domaxis(3) domaxis(4)]);

fontsize=10;

p2=m_contourf(lon_sub,lat_sub,h, 'LevelList',c_ticks,'LineColor', 'none');
hold on;
[C,d]=m_contour(lon_sub,lat_sub,h,contours,'Color','black');
clabel(C,d,contours,'FontSize',8,'Color','black');

%colormap(haxby);
colormap(flip(cmocean('Deep')));
c=colorbar;
ylabel(c,'Depth (m)','FontSize', 14);
shading flat;
cb=colorbar;
ylabel(cb,'Depth (m)');

coastfileplot = 'coastline_f.mat';
m_usercoast(coastfileplot,'patch',[.9 .9 .9]);
m_grid('box','fancy','tickdir','out','fontsize',fontsize);

filename = 'plot_bathy';
print(filename,'-r500','-djpeg')

