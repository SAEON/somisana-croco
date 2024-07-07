%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  make_ERA5.m
% 
%  Create and fill frc and bulk files with ERA5 data.
%  (ERA-5 Reanalysis)
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
%  Updated    D. Donoso, G. Cambon. P. Penven (Oct 2021) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
crocotools_param        
frc_prefix=[frc_prefix,'_ERA5_'];
blk_prefix=[blk_prefix,'_ERA5_'];
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

%
% Get the ERA5 horizontal grids (it should be the same for every month)
%
nc=netcdf([ERA5_dir,'LSM_Y',num2str(Ymin),'M',num2str(Mmin),'.nc']);
disp(['Use this land file :',char([ERA5_dir,'LSM_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'])])

lon1=nc{'lon'}(:);
lat1=nc{'lat'}(:);
[lon1,lat1]=meshgrid(lon1,lat1);

mask=squeeze(nc{'LSM'}(1,:,:)); %we take the first record
mask(mask ~=0 ) = 1;
mask = 1-mask ;
mask(mask ==0 ) = NaN ;
close(nc);


if add_waves == 1
    %
    % get horizontal grid for ocean waves
    % for ERA5 wave variables, resolution change to 0.5° instead of 0.25° 
    % for the atmo variable (0.25° x 0.25° (atmosphere)  vs  0.5° x 0.5° (ocean waves))
    %
    filein=[ERA5_dir,'SWH_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'];
    nc=netcdf(filein);
    lonwave=nc{'lon'}(:);
    latwave=nc{'lat'}(:);
    [lon1wave,lat1wave]=meshgrid(lonwave,latwave);
    maskwave=squeeze(nc{'SWH'}(1,:,:));
    missvalue_wave = ncreadatt(filein,'SWH','missing_value');
    maskwave(maskwave ~= missvalue_wave ) = 1;
    maskwave(maskwave == missvalue_wave ) = NaN ;
    close(nc);
else
    lon1wave = NaN;
    lat1wave = NaN;
    maskwave = NaN;
end


%
%Loop on the years and the months
%
disp(['====================='])
disp(['INTERPOLATION STEP'])
disp(['====================='])
disp(['Loop on the years and the months'])

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
        %-------------------------------------------------------------------%
        %
        % Process time (here in days), with SST file (minimum time step)
        %
        %-------------------------------------------------------------------%
        nc=netcdf([ERA5_dir,'T2M_Y',num2str(Y),'M',num2str(M),'.nc']);
        ERA5_time=nc{'time'}(:);
        close(nc);
        dt=mean(gradient(ERA5_time));
        disp(['dt=',num2str(dt)])
        %-----------------------------------------------------------
        %Variable overlapping timesteps : 2 at the beginning and 2 at the end
        %------------------------------------------------------------
        tlen0=length(ERA5_time);
        disp(['tlen0=',num2str(tlen0)])
        freq=1; % hourly
        itolap=freq*itolap_era5;
        tlen=tlen0+2*itolap;
        disp(['tlen=',num2str(tlen)])
        disp(['Overlap is ',num2str(itolap_era5),' records before and after'])     
        time=0*(1:tlen);
        time(itolap+1:tlen0+itolap)=ERA5_time;   
        disp(['====================='])
        disp('Compute time for croco file')
        disp(['====================='])
        for aa=1:itolap
            time(aa)=time(itolap+1)-(itolap+1-aa)*dt;
        end
        for aa=1:itolap
            time(tlen0+itolap+aa)=time(tlen0+itolap)+aa*dt;
        end
        %-------------------------------------------------------------------%
        %
        % Create the CROCO bulk forcing files
        %
        % ------------------------------------------------------------------%
        %
	disp(['====================='])
        disp('Create the frc/blk netcdf file')
        disp(['====================='])
        blkname=[blk_prefix,'Y',num2str(Y),...
                 'M',num2str(sprintf(Mth_format,M)),nc_suffix];
        frcname=[frc_prefix,'Y',num2str(Y),...
                 'M',num2str(sprintf(Mth_format,M)),nc_suffix];
        if makeblk==1
          disp(['Create a new bulk file: ' blkname])
          create_bulk(blkname,grdname,CROCO_title,time,0);
          nc_add_globatt(blkname,Yorig,Mmin,Dmin,Hmin,Min_min,Smin,'ERA5');
          disp([' '])
        end
        if makefrc==1
          disp(['Create a new forcing file: ' frcname])
          disp([' '])
          create_forcing(frcname,grdname,CROCO_title,...
                          time,0,0,...
                          0,0,0,...
                          0,0,0,0,0,0)
        end
        % %
        % % Add the atmospheric pressure variable to the bulk file
        % %
        if makeblk==1 && add_atm==1
            disp(['Add atmospheric pressure'])
            disp(['========================'])
            nw = netcdf(blkname,'write');
            nw{'patm'}             = ncdouble('bulk_time', 'eta_rho', 'xi_rho');
            nw{'patm'}.long_name   = ncchar('atmospheric pressure above mean seal level');
            nw{'patm'}.long_name   = 'atmospheric pressure above mean seal level';
            nw{'patm'}.units       = ncchar('Pa');
            nw{'patm'}.units       = 'Pa';
            close(nw);
        end
        % % 
        % %
        % % Add the waves
        % %
        if makefrc==1 && add_waves==1
            disp(['Add waves data'])
            disp(['=============='])
        end
	% % 
        % %
        % % Add the tides
        % %
        if makefrc==1 && add_tides==1
             add_tidal_data(tidename,grdname,frcname,Ntides,tidalrank,...
                            Yorig,Y,M,coastfileplot)
        end
        %
        %
        % Open the CROCO forcing files
        if makefrc==1
          nc_frc=netcdf(frcname,'write');
        else
          nc_frc=[];
        end

	% Open the CROCO bulk files
        if makeblk==1
	  nc_blk=netcdf(blkname,'write');
        else
          nc_blk=[];
        end
	%
        % Check if there are ERA5 files for the previous Month
        Mm=M-1;
        Ym=Y;
        if Mm==0
            Mm=12;
            Ym=Y-1;
        end
        %
        fname = [ERA5_dir,'TP_Y',num2str(Ym),'M',num2str(Mm),'.nc'];
        %nc=netcdf([ERA5_dir,'TP_Y',num2str(Ym),'M',num2str(Mm),'.nc']);
        %
        disp(' ')
        disp('======================================================')
        disp('Perform interpolations for the previous month')      
        disp('======================================================')
        disp(' ')
        if exist(fname)==0
            disp(['No data for the previous month: using current month'])
            tndx=1;
            Mm=M;
            Ym=Y;
        else
            nc=netcdf(fname);
            tndx=length(nc('time'));
            if makefrc==1
               for aa=1:itolap
                 nc_frc{'sms_time'}(aa)=nc{'time'}(tndx-(itolap-aa));
               end
            end
            %
	    if makeblk==1
               for aa=1:itolap
                 nc_blk{'bulk_time'}(aa)=nc{'time'}(tndx-(itolap-aa));
               end
	    end
            close(nc)
        end
        %
        % Perform interpolations for the previous month or repeat the first one
        %
        for aa=1:itolap
            if exist(fname)==0
                aa0=aa;
            else
                aa0=tndx-(itolap-aa);
            end
            interp_ERA5(ERA5_dir,Ym,Mm,Roa,interp_method,lon1,lat1,lon1wave,lat1wave, ...
                        mask,maskwave,aa0,nc_frc,nc_blk,lon,lat,angle,aa,add_waves,add_atm)
        end  
        %######################################################################      
        %   
        disp(' ')
        disp('======================================================')
        disp('Perform interpolations for the current month')
        disp('======================================================')
        disp(' ')
        
        % Perform interpolations for the current month
        %
        
        for tndx=1:tlen0
            if mod(tndx,6)==0
                disp(['Step: ',num2str(tndx),' of ',num2str(tlen0)])
            end
            interp_ERA5(ERA5_dir,Y,M,Roa,interp_method,lon1,lat1,lon1wave,lat1wave,...
                        mask,maskwave,tndx,nc_frc,nc_blk,lon,lat,angle,tndx+itolap,add_waves,add_atm)	
        end
        
        disp(' ')      
        disp('======================================================')
        disp('Perform interpolations for next month')    
        disp('======================================================')
        disp(' ')
        %######################################################################
        % Read ERA5 file for the next month
        %
        Mp=M+1;
        Yp=Y;
        if Mp==13
            Mp=1;
            Yp=Y+1;
        end
        
        fname=[ERA5_dir,'TP_Y',num2str(Yp),'M',num2str(Mp),'.nc'];
        
        if exist(fname)==0
            disp(['No data for the next month: using current month'])
            tndx=tlen0;
            Mp=M;
            Yp=Y;
        else
            nc=netcdf(fname);

            if makefrc==1
              disp('sms_time')
              for tndx=tlen0+itolap+1:tlen;
                nc_frc{'sms_time'}(tndx)=nc{'time'}(tndx-tlen0-(itolap));
              end;
            end

	    if makeblk==1
              disp('bulk_time')
              for tndx=tlen0+itolap+1:tlen;
                nc_blk{'bulk_time'}(tndx)=nc{'time'}(tndx-tlen0-(itolap));
              end
	    end
            close(nc)
        end
        %
        % Perform the interpolations for the next month
        %
        disp('Last steps')
        for tndx=tlen0+itolap+1:tlen;
            disp(['tndx= ',num2str(tndx)])
            tout=tndx;
            disp(['tout=tndx ',num2str(tndx)])
            if Mp==M
                %tin=tlen0;       % persistency if current month is used
                tin=tndx-2*itolap
                disp(['tin=',num2str(tin)])
            else
                tin=tndx-tlen0-itolap;
                disp(['tin=',num2str(tin)])
            end
            interp_ERA5(ERA5_dir,Yp,Mp,Roa,interp_method,lon1,lat1,lon1wave,lat1wave,...
                        mask,maskwave,tin,nc_frc,nc_blk,lon,lat,angle,tout,add_waves,add_atm)           
        end;
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
%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%
disp('======================================================')
if SPIN_Long>0
disp('Add spin up phase')
    M=Mmin-1;
    Y=Ymin-SPIN_Long;
    for month=1:12*SPIN_Long
        M=M+1;
        if M==13
            M=1; 
            Y=Y+1;
        end
        %
        % Bulk files
        %
        %
        % Copy the file
        %
        blkname=[blk_prefix,'Y',num2str(Ymin),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
        blkname2=[blk_prefix,'Y',num2str(Y),'M',num2str(sprintf(Mth_format,M)),nc_suffix];
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
%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
    disp(' ')
    disp('======================================================')
    disp(' Make a few plots...')
    slides=[1 12 24 36];
    
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

