%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  make_CSAG.m
% 
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%  script is based on make_ECMWF.m and make_CFSR.m as part of croco_tools, 
%  adapted to process data provided by the Climate Systems Analysis Group
%  (CSAG) at UCT
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
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
% 
%  Adapted from a previous verions from
%  Alvaro Peliz (U. Aveiro) & Patrick Marchesiello (IRD) - 2005
%  Created by Serena Illig 2010
%
%  Created    January 2016 (S. Illig and E. Cadier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param
%
  frc_prefix=[frc_prefix,'_CSAG_'];                           
  blk_prefix=[blk_prefix,'_CSAG_'];
%
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
disp(' ')
disp([' Read the grid in ',grdname])
nc=netcdf(grdname);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
lonu=nc{'lon_u'}(:);
latu=nc{'lat_u'}(:);
lonv=nc{'lon_v'}(:);
latv=nc{'lat_v'}(:);
angle=nc{'angle'}(:);
close(nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data from netcdf format saved on my computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if Download_data==1
%
% This step is skipped as we already extract a subset of the CSAG dataset 
% using cdo when downloading from /terra/data/windatlas/wasa/netcdf/3km
%
% end

if makefrc==1 | makeblk==1
  %
  % Get the CSAG horizontal grid
  %
  fname=[My_CSAG_dir,'/wrfout_d03_wind_10m_regridded_sub.nc'];
  nc=netcdf(fname); % any file will do - grid same for all
  lon1=squeeze(nc{'lon'}(:));
  lat1=squeeze(nc{'lat'}(:));
  [lon1,lat1]=meshgrid(lon1,lat1);
  close(nc);  
  %
  % Get the CSAG time array
  %
  disp('==========================')
  disp('getting generic time array...')
  disp('==========================')
  %
  % this is the same for all variables except 'rain' (which has a few steps
  % missing at the beginning). So we need to get two time arrays- one for 
  % rain and one for all the other variables
  %
  time=ncread(fname,'time');
  time=time/24; % days since reference time
  ref_time=datenum(1989,1,1);
  time=time+ref_time; % time as matlab datenum
  datestrings=datestr(time);
  datetimes=datetime(datestrings);
  time=time-datenum(Yorig,1,1); % time is now in days since 1-Jan-Yorig
  %
  disp('==========================')
  disp('getting rain time array...')
  disp('==========================')
  %
  time_rain=ncread([My_CSAG_dir,'/wrfout_d03_rain_regridded_sub.nc'],'time');
  time_rain=time_rain/24; % days since reference time
  time_rain=time_rain+ref_time; % time as matlab datenum
  datestrings_rain=datestr(time_rain);
  datetimes_rain=datetime(datestrings_rain);
  time_rain=time_rain-datenum(Yorig,1,1); % time is now in days since 1-Jan-Yorig
  %
  % Loop on the years and the months
  %
   for Y=Ymin:Ymax
    if Y==Ymin 
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
      mo_max=12;
    end
     for M=mo_min:mo_max
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
            ' - month ',num2str(M)])
      disp(' ')
      %
      % Get the time subset for this month
      %
      % specify start time
      datenum_start=datenum(Y,M,1,0,0,0)-itolap_csag;
      datestr_start=datestr(datenum_start,0);
      datetime_start=datetime(datestr_start);
      % specify end time
      nmax=daysinmonth(Y,M);
      datenum_end=datenum(Y,M,nmax,23,0,0)+ itolap_csag;
      datestr_end=datestr(datenum_end,0);
      datetime_end=datetime(datestr_end);
      % get the time indices for this month
      indx_var=find(isbetween(datetimes,datetime_start,datetime_end));
      time_sub=time(indx_var);
      % noticed 1 day of duplicate time-steps so remove these (this is
      % right in the beginning of the file - starting tstep 240 - so we
      % could just avoid using this time period rather than coding around
      % it... oh well, no harm in including the next couple of lines)
      [time_sub, indx_sub, ~] = unique(time_sub);
      indx_var=indx_var(indx_sub); % this gets the indices used to extract the variable from the original data
      %
      % rain time array will have different indices than the other
      % variables
      indx_rain=find(isbetween(datetimes_rain,datetime_start,datetime_end));
      % no need to remove duplicates for the rain time array (checked)
      % 
      %
      %-------------------------------------------------------------------%
      %
      % Create the CROCO forcing files
      %
      % ------------------------------------------------------------------%
      %
      disp(['====================='])
      disp('Create the frc/blk netcdf file')
      disp(['====================='])
      %
      blkname=[blk_prefix,'Y',num2str(Y),...
               'M',num2str(M),nc_suffix];
      frcname=[frc_prefix,'Y',num2str(Y),...
               'M',num2str(M),nc_suffix];
      if makeblk==1
        disp(['Create a new bulk file: ' blkname])
        create_bulk(blkname,grdname,CROCO_title,time_sub,0);
	    disp([' '])
      end
      if makefrc==1
        disp(['Create a new forcing file: ' frcname])
	    disp([' '])
        create_forcing(frcname,grdname,CROCO_title,...
                       time_sub,0,0,...
                       0,0,0,...
                       0,0,0,0,0,0)
      end
      %
      % Open the CROCO forcing files
      if makefrc==1
        nc_frc=netcdf(frcname,'write');
      else
        nc_frc=[];
      end
      if makeblk==1
        nc_blk=netcdf(blkname,'write');
      else
        nc_blk=[];
      end
      %
	  
      %######################################################################      
      %   
      disp(' ')
      disp('======================================================')
      disp('Perform interpolations for the current month')
      disp('======================================================')
      disp(' ')

      % Perform interpolations for the current month
      %

      for tndx=1:length(time_sub) %
        if mod(tndx,20)==0
          disp(['Step: ',num2str(tndx),' of ',num2str(length(time_sub))])
        end
	    interp_CSAG(My_CSAG_dir,Roa,interp_method,lon1,lat1,...
		    tndx,indx_var,indx_rain,nc_frc,nc_blk,lon,lat,angle)	
      end
      
      %
      % Close the CROCO forcing files
      %
      if ~isempty(nc_frc)
        close(nc_frc);
      end
      if ~isempty(nc_blk)
        close(nc_blk);
      end
    end
   end
 end
%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%
disp('======================================================')
disp('Add spin up phase')      
if SPIN_Long>0
  M=Mmin-1;
  Y=Ymin-SPIN_Long;
  for month=1:12*SPIN_Long
    M=M+1;
    if M==13
      M=1; 
      Y=Y+1;
    end
    %
    % Forcing files
    %
    if makefrc==1
      %
      % Copy the file
      %
      frcname=[frc_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      frcname2=[frc_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',frcname2]) 
      eval(['!cp ',frcname,' ',frcname2]) 
      %
      % Change the time
      %
      nc=netcdf(frcname2,'write');
      time=nc{'sms_time'}(:)-365.*(Ymin-Y);%+datenum(Yorig,1,1);
      %[y,m,d,h,mi,s]=datevec(time);
      %dy=Ymin-Y;
      %y=y-dy;
      %time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
      nc{'sms_time'}(:)=time;
      close(nc)
    end
    %
    % Bulk files
    %
    if makeblk==1
      %
      % Copy the file
      %
      blkname=[blk_prefix,'Y',num2str(Ymin),'M',num2str(M),nc_suffix];
      blkname2=[blk_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
      disp(['Create ',blkname2]) 
      eval(['!cp ',blkname,' ',blkname2]) 
      %
      % Change the time
      %
      nc=netcdf(blkname2,'write');
      time=nc{'bulk_time'}(:)-365.*(Ymin-Y);%+datenum(Yorig,1,1);
      %[y,m,d,h,mi,s]=datevec(time);
      %dy=Ymin-Y;
      %y=y-dy;
      %time=datenum(y,m,d,h,mi,s)-datenum(Yorig,1,1);
      nc{'bulk_time'}(:)=time;
      close(nc)
    end
  end
end
%---------------------------------------------------------------
% Create MC file
%---------------------------------------------------------------
%if createmc==1
%    if makefrc
%       make_mc_forcing('CFSR',21600.,itolap_ncep);
%   end
%   if makeblk
%       make_mc_bulk('CFSR',21600.,itolap_ncep);
%   end
%end
%
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp('======================================================')
  disp(' Make a few plots...')
  slides=[1 12 24 36];
  if makeblk
    test_forcing(blkname,grdname,'tair',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'rhum',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'prate',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'uwnd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'vwnd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'wspd',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'radlw',slides,3,coastfileplot)
    figure
    test_forcing(blkname,grdname,'radsw',slides,3,coastfileplot)
  end
  if makefrc
    test_forcing(frcname,grdname,'sustr',slides,3,coastfileplot)
    figure
    test_forcing(frcname,grdname,'svstr',slides,3,coastfileplot)  
  end
end

