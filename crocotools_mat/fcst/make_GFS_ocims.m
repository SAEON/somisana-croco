function make_GFS_ocims(NY,NM,ND,delta_days,hdays,fdays)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% G Fearon Oct 2020:
% This script is adapted from make_GFS.m as part of croco_tools, but here we
% make it a function to be called from within our over-arching python
% script, and we don't use matlab's opendap to download the GFS data - that
% gets done in a separate previous step. So instead of calling
% download_GFS.m we call our own script extract_GFS.m which gets data from
% already downloaded grib files
%
% Create and fill frc and bulk files with GFS data.
% for a forecast run
%
% The on-line reference to GFS is at
% http://nomad3.ncep.noaa.gov/
% 
%  Further Information:  
%  http://www.croco-ocean.org
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
%  Copyright (c) 2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    9-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Updated    12-Feb-2016 by P. Marchesiello
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param;

%gfs_grb_dir=[DATADIR,'gfs/']; % Changed to DATADIR all downloaded enviromental data now stored in the same place
%
frc_prefix=[frc_prefix,'_GFS_'];
blk_prefix=[blk_prefix,'_GFS_'];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% time (in matlab time)
%
today=datenum(NY,NM,ND,0,0,0);
%
% GFS data name
%
gfs_name=[FORC_DATA_DIR,'GFS_',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),'.nc'];
%
% Extract data from the already downloaded enviromental data files and write 
% into a more croco_tools friendly format. Store the output in the forcing temp
% directory created everyday.
extract_GFS(DATADIR,today,delta_days,hdays,fdays,Yorig,gfs_name)
%
if level==0
  nc_suffix='.nc';
else
  nc_suffix=['.nc.',num2str(level)];
  grdname=[grdname,'.',num2str(level)];
end
%
% Get the model grid
%
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)
cosa=cos(angle);
sina=sin(angle);
%
% Get the GFS grid 
% 
nc=netcdf(gfs_name);
lon1=nc{'lon'}(:);
lat1=nc{'lat'}(:);
time=nc{'time'}(:);
mask=nc{'mask'}(:);
tlen=length(time);
%
% bulk and forcing files
%
blkname=[blk_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),nc_suffix];
disp(['Create a new bulk file: ' blkname])
create_bulk(blkname,grdname,CROCO_title,time,0);
nc_blk=netcdf(blkname,'write');
%
% Commenting out the making and filling of the forcing file, as our
% forecasts use bulk surface forcing for all atmospheric variables
%
% frcname=[frc_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),nc_suffix];
% disp(['Create a new forcing file: ' frcname])
% create_forcing(frcname,grdname,CROCO_title,...
%                        time,0,0,...
%                        0,0,0,...
%   	               0,0,0,0,0,0)
% nc_frc=netcdf(frcname,'write');
%
% Loop on time
%
missval=nan;
default=nan;
for l=1:tlen
  %disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
  var=squeeze(nc{'tair'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'tair'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'tair'}(l-1,:,:)); 
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'tair'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end

  var=squeeze(nc{'rhum'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'rhum'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'rhum'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'rhum'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
  
  var=squeeze(nc{'prate'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'prate'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'prate'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'prate'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
  
  var=squeeze(nc{'wspd'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'wspd'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'wspd'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'wspd'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
  
  %Zonal wind speed
  var=squeeze(nc{'uwnd'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    uwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
  else
    var=squeeze(nc{'uwnd'}(l-1,:,:));
    uwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
  end
  
  %Meridian wind speed
  var=squeeze(nc{'vwnd'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    vwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
  else
    var=squeeze(nc{'vwnd'}(l-1,:,:));
    vwnd=get_missing_val(lon1,lat1,var,missval,Roa,default);
    vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
  end
  
%   nc_frc{'uwnd'}(l,:,:)=rho2u_2d(uwnd.*cosa+vwnd.*sina);
%   nc_frc{'vwnd'}(l,:,:)=rho2v_2d(vwnd.*cosa-uwnd.*sina);
  
  nc_blk{'uwnd'}(l,:,:)=rho2u_2d(uwnd.*cosa+vwnd.*sina);
  nc_blk{'vwnd'}(l,:,:)=rho2v_2d(vwnd.*cosa-uwnd.*sina);
  
  %Net longwave flux
  var=squeeze(nc{'radlw'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radlw'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'radlw'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radlw'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
  
  %Downward longwave flux
  var=squeeze(nc{'radlw_in'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radlw_in'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'radlw_in'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radlw_in'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
   
  %Net solar short wave radiation
  var=squeeze(nc{'radsw'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radsw'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  else
    var=squeeze(nc{'radsw'}(l-1,:,:));
    var=get_missing_val(lon1,lat1,var,missval,Roa,default);
    nc_blk{'radsw'}(l,:,:)=interp2(lon1,lat1,var,lon,lat,interp_method);
  end
  
  var=squeeze(nc{'tx'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    tx=get_missing_val(lon1,lat1,var,missval,Roa,default);
    tx=interp2(lon1,lat1,tx,lon,lat,interp_method);
  else
    var=squeeze(nc{'tx'}(l-1,:,:));
    tx=get_missing_val(lon1,lat1,var,missval,Roa,default);
    tx=interp2(lon1,lat1,tx,lon,lat,interp_method);
  end
  
  var=squeeze(nc{'ty'}(l,:,:));
  if mean(mean(isnan(var)~=1))
    ty=get_missing_val(lon1,lat1,var,missval,Roa,default);
    ty=interp2(lon1,lat1,ty,lon,lat,interp_method);
  else
    var=squeeze(nc{'ty'}(l-1,:,:));
    ty=get_missing_val(lon1,lat1,var,missval,Roa,default);
    ty=interp2(lon1,lat1,ty,lon,lat,interp_method);
  end
   
%   nc_frc{'sustr'}(l,:,:)=rho2u_2d(tx.*cosa+ty.*sina);
%   nc_frc{'svstr'}(l,:,:)=rho2v_2d(ty.*cosa-tx.*sina);
  
  nc_blk{'sustr'}(l,:,:)=rho2u_2d(tx.*cosa+ty.*sina);
  nc_blk{'svstr'}(l,:,:)=rho2v_2d(ty.*cosa-tx.*sina);
end
% 
%close(nc_frc);
close(nc_blk);
close(nc)
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
% if makeplot==1
%   disp(' ')
%   disp(' Make a few plots...')
%   slides=[10 12 14 16]; 
%   test_forcing(blkname,grdname,'tair',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'rhum',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'prate',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'wspd',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'radlw',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'radlw_in',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'sustr',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'svstr',slides,3,coastfileplot)
%   figure
%   test_forcing(blkname,grdname,'radsw',slides,3,coastfileplot)
% end
end








