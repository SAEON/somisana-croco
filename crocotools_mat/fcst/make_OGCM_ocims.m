function make_OGCM_ocims(NY,NM,ND,hdays,makeini)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% G Fearon Oct 2020:
% This script is adapted from make_OGCM_frcst.m as part of croco_tools, but here we
% make it a function to be called from within our over-arching python
% script. The mercator data is downloaded in a separate previous step.
%
% Create and fill CROCO clim and bry files with OGCM data.
% for a forecast run
%
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
%  Updated    8-Sep-2006 by Pierrick Penven
%  Updated    20-Aug-2008 by Matthieu Caillaud & P. Marchesiello
%  Updated   12-Feb-2016 by P. Marchesiello
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
%close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param;
%
% name of raw downloaded enviromental data file for this day. This file should be in DATADIR 
% where all downloaded enviromental data is stored.
MERCATOR_name_raw=[DATADIR,'mercator_',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),'.nc'];
%
% write it into a more croco_tools friendly format. The output goes to forcing file folder FORC_DATA_DIR
% where dynamic forcing files are created
MERCATOR_name=[FORC_DATA_DIR,OGCM_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),'.cdf'];
write_mercator_ocims(MERCATOR_name,MERCATOR_name_raw,Yorig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of user input  parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
% Get the OGCM grid 
% 
nc=netcdf(MERCATOR_name);
lonT=nc{'lonT'}(:);
latT=nc{'latT'}(:);
lonU=nc{'lonU'}(:);
latU=nc{'latU'}(:);
lonV=nc{'lonV'}(:);
latV=nc{'latV'}(:);
Z=-nc{'depth'}(:);
NZ=length(Z);
NZ=NZ-rmdepth;
Z=Z(1:NZ);
%
% make_OGCM_frcst.m subsets the time to only three time steps... 
% not sure why? Probably just to save time... Our model is small and we've 
% set up the download so that the file will contain data which more than 
% covers the simulation period. So I'm just going to write every timestep 
% as saved in the raw OGCM file
time=nc{'time'}(:); % this time is days since datenum(Yorig,1,1)
time_cycle=0;

close(nc)

%
% Initial file 
%
if makeini==1
  %
  % by definition our simulation starts at NY-NM-ND 00:00:00 - hdays
  tini=datenum(NY,NM,ND)-hdays-datenum(Yorig,1,1); 
  %
  % Mercator data is saved at time 12:00:00 daily, so I should do some 
  % temporal interpolation to get the interpolated Mercator data at tini.
  % Beind lazy here, I'm just using the nearest indx to tini. 
  % This might cause stability issues at the first initialisation, but
  % thereafter we're initialising from the previous model solution, so I'm
  % guessing this is no big deal
  [~,tndx]=min(abs(time-tini)); % using the nearest indx to the start time
  ininame=[ini_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),nc_suffix];
  %disp(['Create an initial file for ',num2str(NY),num2str(NM),num2str(ND);])
  create_inifile(ininame,grdname,CROCO_title,...
                 theta_s,theta_b,hc,N,...
                 tini,'clobber',vtransform);
  nc_ini=netcdf(ininame,'write');
  
  interp_OGCM_ocims(MERCATOR_name,Roa,interp_method,...
              lonU,latU,lonV,latV,lonT,latT,Z,tndx,...
              nc_ini,[],lon,lat,angle,h,1,vtransform)
  close(nc_ini)
end

% Clim and Bry files
%
if makeclim==1 | makebry==1
  if makebry==1
    bryname=[bry_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),nc_suffix];
    create_bryfile(bryname,grdname,CROCO_title,[1 1 1 1],...
                   theta_s,theta_b,hc,N,...
                   time,time_cycle,'clobber',vtransform);
    nc_bry=netcdf(bryname,'write');
  else
    nc_bry=[];
  end
  if makeclim==1
    clmname=[clm_prefix,num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),nc_suffix];
    create_climfile(clmname,grdname,CROCO_title,...
                    theta_s,theta_b,hc,N,...
                    time,time_cycle,'clobber',vtransform);
    nc_clm=netcdf(clmname,'write');
  else
    nc_clm=[];
  end

%
% Perform the interpolations for all selected records
%
for tndx=1:length(time)
  %disp([' Time step : ',num2str(tndx),' of ',num2str(length(time)),' :'])
  interp_OGCM_ocims(MERCATOR_name,Roa,interp_method,...
                    lonU,latU,lonV,latV,lonT,latT,Z,tndx,...
		    nc_clm,nc_bry,lon,lat,angle,h,tndx,vtransform)
end

%
% Close the CROCO files
%
  if ~isempty(nc_clm)
    close(nc_clm);
  end
  if ~isempty(nc_bry)
    close(nc_bry);
  end
%
end
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
% if makeplot==1
%   disp(' ')
%   disp(' Make a few plots...')
%   if makeini==1
%     ininame=[ini_prefix,num2str(rundate),nc_suffix];
%     figure
%     test_clim(ininame,grdname,'temp',1,coastfileplot)
%     figure
%     test_clim(ininame,grdname,'salt',1,coastfileplot)
%   end
%   if makeclim==1
%     clmname=[clm_prefix,num2str(rundate),nc_suffix];
%     figure
%     test_clim(clmname,grdname,'temp',1,coastfileplot)
%     figure
%     test_clim(clmname,grdname,'salt',1,coastfileplot)
%   end
%   if makebry==1
%     bryname=[bry_prefix,num2str(rundate),nc_suffix];
%     figure
%     test_bry(bryname,grdname,'temp',1,obc)
%     figure
%     test_bry(bryname,grdname,'salt',1,obc)
%     figure
%     test_bry(bryname,grdname,'u',1,obc)
%     figure
%     test_bry(bryname,grdname,'v',1,obc)
%   end
% end

end
