function reformat_GFS(Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  G Fearon Mar 2024:
%  Reformat the GFS data in grb format into a format which can be used to make blk croco input files
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  %
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the data used in the download of the GFS data
% assuming this script is run in the same directory where the data were downloaded 
gfs_env=[pwd,'/gfs.env']
fileID = fopen(gfs_env, 'r');
[run_date,delta_days,hdays,fdays]=read_gfs_env(gfs_env);
%
% Extract year, month, day, hour
[NY, NM, ND, NH] = datevec(run_date);

% filename for reformatted GFS data
%
gfs_name=[pwd,'/for_croco/','GFS_',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),'_',num2str(NH,'%02.f'),'.nc'];

% set up toolbox for reading grib files
setup_nctoolbox

% checked that these start and end times are correct, i.e. tested against the
% script reads the data. There is also a check in the file which will alert
% us to a problem with this
%
% looking 6 hours either side of the simulation window to make sure our
% forcing completely covers the simulation time (this is mirrored in the 
% download of the data)
% I'm now commenting this as we are now already reading in the extended values used in the GFS download
%hdays=hdays+0.25;
%fdays=fdays+0.25;

datenum_start=datenum(run_date-hdays);
datenum_end=datenum(run_date+fdays);

gfstime=datenum_start+1/24:1/24:datenum_end; % start from datenum_start+1/24 as the first time step is the first forecast hour from datenum_start
gfstime_check=gfstime; % gfstime is overwritten in this script and checked at the end to make sure our times were correctly assigned

% use any file to get the gfs grid
dv=datevec(datenum_start);
NY=dv(1);NM=dv(2);ND=dv(3);NH=dv(4);
%NY=year(datenum_start);NM=month(datenum_start);ND=day(datenum_start);NH=hour(datenum_start);
gfs_grb_name=[pwd,'/',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),num2str(NH,'%02.f'),'_f001.grb'];
nco=ncgeodataset(gfs_grb_name);
lat=nco{'lat'}(:);
lon=nco{'lon'}(:);
%nco.variables
mask=squeeze(nco{'Land_cover_0__sea_1__land_surface'}(:));
mask(mask==1)=NaN;
mask(isfinite(mask))=1;
clear nco;

% Initialisation
N=length(gfstime);
[M,L]=size(mask);
tx=zeros(N,M,L);
ty=tx;
tair=tx;
rhum=tx;
prate=tx;
wspd=tx;
radlw=tx;
radsw=tx;
radlw_in=tx;
uwnd=tx;
vwnd=tx;

%
%==================================================
% Get the data for the hindcast portion of the run
%==================================================
%
n=0;
% this looping method mirrors that used to download the data
datenum_hist=datenum_start;
datenum_latest=run_date+delta_days;
while datenum_hist<datenum_latest
    %disp(['reading historical GFS data for ',datestr(datenum_hist)]);
    dv=datevec(datenum_hist);
    NY=dv(1);NM=dv(2);ND=dv(3);NH=dv(4);
    %NY=year(datenum_hist);NM=month(datenum_hist);ND=day(datenum_hist);NH=hour(datenum_hist);
    for frcst=1:6
        n=n+1;
        gfstime(n)=datenum(NY,NM,ND,NH,0,0)+frcst/24;
        gfs_grb_name=[pwd,'/',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),num2str(NH,'%02.f'),'_f',num2str(frcst,'%03.f'),'.grb'];
        
        nco=ncgeodataset(gfs_grb_name);
        %
        %nco.variables
        %nco.variable(['Momentum_flux_u-component_surface_',num2str(frcst),'_Hour_Average']).attributes
        %
        % read each variable and populate the output arrays
        %
        u_n=squeeze(nco{'u-component_of_wind_height_above_ground'}(:));
        v_n=squeeze(nco{'v-component_of_wind_height_above_ground'}(:));
        tair_n=squeeze(nco{'Temperature_height_above_ground'}(:));
        rhum_n=squeeze(nco{'Relative_humidity_height_above_ground'}(:));
        prate_n=squeeze(nco{'Precipitation_rate_surface'}(:));
        % A bunch of variables are (rather annoyingly) written out as the 
        % accumulated average over each 6 hour forecast period: 
        tx_n_ave=squeeze(nco{['Momentum_flux_u-component_surface_',num2str(frcst),'_Hour_Average']}(:));
        ty_n_ave=squeeze(nco{['Momentum_flux_v-component_surface_',num2str(frcst),'_Hour_Average']}(:));
        dradlw_n_ave=squeeze(nco{['Downward_Long-Wave_Radp_Flux_surface_',num2str(frcst),'_Hour_Average']}(:));
        uradlw_n_ave=squeeze(nco{['Upward_Long-Wave_Radp_Flux_surface_',num2str(frcst),'_Hour_Average']}(:));
        dradsw_n_ave=squeeze(nco{['Downward_Short-Wave_Radiation_Flux_surface_',num2str(frcst),'_Hour_Average']}(:));
        uradsw_n_ave=squeeze(nco{['Upward_Short-Wave_Radiation_Flux_surface_',num2str(frcst),'_Hour_Average']}(:));
        % but we want to convert these accumulated averages into individual one-hour averages 
        % see FAQ "How can the individual one-hour averages be computed" from https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html
        % Excerpt from there:
        % You can compute the one-hour average (X) ending at hour N by using the N-hour average (a) and the (N-1)-hour average (b) as follows:
        % X = N*a - (N-1)*b
        % So if you want the 1-hour Average for the period initial+3 to initial+4 (X), you would use the 4-hour Average (initial+0 to initial+4) as (a) and the 3-hour Average (initial+0 to initial+3) as (b) as follows:
        % X = 4*a - 3*b
        if frcst==1
            tx_n=tx_n_ave;
            ty_n=ty_n_ave;
            dradlw_n=dradlw_n_ave;
            uradlw_n=uradlw_n_ave;
            dradsw_n=dradsw_n_ave;
            uradsw_n=uradsw_n_ave;
        else
            tx_n=tx_n_ave*(frcst)-tx_n_ave_prev*(frcst-1);
            ty_n=ty_n_ave*(frcst)-ty_n_ave_prev*(frcst-1);
            dradlw_n=dradlw_n_ave*(frcst)-dradlw_n_ave_prev*(frcst-1);
            uradlw_n=uradlw_n_ave*(frcst)-uradlw_n_ave_prev*(frcst-1);
            dradsw_n=dradsw_n_ave*(frcst)-dradsw_n_ave_prev*(frcst-1);
            uradsw_n=uradsw_n_ave*(frcst)-uradsw_n_ave_prev*(frcst-1);
        end
        
        % save accumulated average data for use in the next time-step
        tx_n_ave_prev=tx_n_ave;
        ty_n_ave_prev=ty_n_ave;
        dradlw_n_ave_prev=dradlw_n_ave;
        uradlw_n_ave_prev=uradlw_n_ave;
        dradsw_n_ave_prev=dradsw_n_ave;
        uradsw_n_ave_prev=uradsw_n_ave;
        
        % Convert the units and write to output variables
        %
        %1: Air temperature: Convert from Kelvin to Celsius
        tair(n,:,:)=tair_n-273.15;
        %
        % 2: Relative humidity: Convert from % to fraction
        rhum(n,:,:)=rhum_n/100;
        %
        % 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
        prate_n=prate_n*0.1*(24*60*60.0);
        prate_n(abs(prate_n)<1.e-4)=0;
        prate(n,:,:)=prate_n;
        %
        % 4: Net shortwave flux: [W/m^2]
        %    CROCO convention: positive downward: same as GFS
        % ?? albedo ??
        radsw_n=dradsw_n - uradsw_n;
        radsw_n(radsw_n<1.e-10)=0;
        radsw(n,:,:)=radsw_n;
        %
        % 5: Net outgoing Longwave flux:  [W/m^2]
        %    CROCO convention: positive upward (opposite to nswrs)
        %    GFS convention: positive downward --> * (-1)
        %    input: downward longwave rad. and
        %    skin temperature.
        radlw(n,:,:) = uradlw_n - dradlw_n;
        radlw_in(n,:,:)=dradlw_n;
        %
        % 6: Wind speed
        wspd(n,:,:)=sqrt(u_n.^2+v_n.^2);
        %
        % 7:  Wind vectors
        uwnd(n,:,:)=u_n;
        vwnd(n,:,:)=v_n;
        % 
        % 7: Wind stress
        tx(n,:,:)=tx_n;
        ty(n,:,:)=ty_n;
        
        clear nco;
    end
    datenum_hist=datenum_hist+6/24;
end

%
%==================================================
% Get the data for the forecast portion of the run
%==================================================
%
% A lot of the code for extracting the data over the forecast portion
% is repeated from the hindcast portion, so I really should have put it
% into another function, maybe called get_GFS_ocims.m, but I got lazy.
% Also, you have to have data from the previous time-step to process
% some of the variables, which complicates a straight extraction of data 
% for a specific time-step in a standalone function. 
% Anyway, it all seems to work so leaving it for now...
%
disp(['reading forecast GFS data from ',datestr(datenum_latest)]);

dv=datevec(datenum_latest);
NY=dv(1);NM=dv(2);ND=dv(3);NH=dv(4);
%NY=year(datenum_latest);NM=month(datenum_latest);ND=day(datenum_latest);NH=hour(datenum_latest);
fhours=(fdays-delta_days)*24;
for frcst=1:fhours
    n=n+1;
    gfstime(n)=datenum(NY,NM,ND,NH,0,0)+frcst/24;
    gfs_grb_name=[pwd,'/',num2str(NY),num2str(NM,'%02.f'),num2str(ND,'%02.f'),num2str(NH,'%02.f'),'_f',num2str(frcst,'%03.f'),'.grb'];
    
    % we need to check if the grb file exists here as data is only stored 
    % in 3 hour increments from frcst=120 hrs onwards. In this case we just 
    % use data from the last available time i.e. we skip the part
    % where we read data from the file
    if exist(gfs_grb_name, 'file')
                
        nco=ncgeodataset(gfs_grb_name);
        %nco.variables
        %
        % read each variable and populate the output arrays
        %
        u_n=squeeze(nco{'u-component_of_wind_height_above_ground'}(:));
        v_n=squeeze(nco{'v-component_of_wind_height_above_ground'}(:));
        tair_n=squeeze(nco{'Temperature_height_above_ground'}(:));
        rhum_n=squeeze(nco{'Relative_humidity_height_above_ground'}(:));
        prate_n=squeeze(nco{'Precipitation_rate_surface'}(:));
        % A bunch of variables are (rather annoyingly) written out as the
        % accumulated average over each 6 hour forecast period:
        % averaging period is reset every 6 hours of the forecast, so we didn't
        % need the next line in reading the historical data, but as we go past
        % 6 hourse in the forecast we need to do this:
        frcst_ave=mod(frcst-1,6)+1; % (tested this works, otherwise next lines result in an error)
        tx_n_ave=squeeze(nco{['Momentum_flux_u-component_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        ty_n_ave=squeeze(nco{['Momentum_flux_v-component_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        dradlw_n_ave=squeeze(nco{['Downward_Long-Wave_Radp_Flux_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        uradlw_n_ave=squeeze(nco{['Upward_Long-Wave_Radp_Flux_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        dradsw_n_ave=squeeze(nco{['Downward_Short-Wave_Radiation_Flux_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        uradsw_n_ave=squeeze(nco{['Upward_Short-Wave_Radiation_Flux_surface_',num2str(frcst_ave),'_Hour_Average']}(:));
        % but we want to convert these accumulated averages into individual one-hour averages
        % see FAQ "How can the individual one-hour averages be computed" from https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html
        % Excerpt from there:
        % You can compute the one-hour average (X) ending at hour N by using the N-hour average (a) and the (N-1)-hour average (b) as follows:
        % X = N*a - (N-1)*b
        % So if you want the 1-hour Average for the period initial+3 to initial+4 (X), you would use the 4-hour Average (initial+0 to initial+4) as (a) and the 3-hour Average (initial+0 to initial+3) as (b) as follows:
        % X = 4*a - 3*b
        if (frcst<=120 && frcst_ave==1) || (frcst>120 && frcst_ave == 3) % here we're handling data going to 3 hourly after frcst = 120 hrs
            tx_n=tx_n_ave;
            ty_n=ty_n_ave;
            dradlw_n=dradlw_n_ave;
            uradlw_n=uradlw_n_ave;
            dradsw_n=dradsw_n_ave;
            uradsw_n=uradsw_n_ave;
        else
            % note here we have to use frcst_ave, not frcst 
            tx_n=tx_n_ave*(frcst_ave)-tx_n_ave_prev*(frcst_ave-1);
            ty_n=ty_n_ave*(frcst_ave)-ty_n_ave_prev*(frcst_ave-1);
            dradlw_n=dradlw_n_ave*(frcst_ave)-dradlw_n_ave_prev*(frcst_ave-1);
            uradlw_n=uradlw_n_ave*(frcst_ave)-uradlw_n_ave_prev*(frcst_ave-1);
            dradsw_n=dradsw_n_ave*(frcst_ave)-dradsw_n_ave_prev*(frcst_ave-1);
            uradsw_n=uradsw_n_ave*(frcst_ave)-uradsw_n_ave_prev*(frcst_ave-1);
        end
        % save accumulated average data for use in the next time-step
        tx_n_ave_prev=tx_n_ave;
        ty_n_ave_prev=ty_n_ave;
        dradlw_n_ave_prev=dradlw_n_ave;
        uradlw_n_ave_prev=uradlw_n_ave;
        dradsw_n_ave_prev=dradsw_n_ave;
        uradsw_n_ave_prev=uradsw_n_ave;
        
        % Transform the variables and write to output variables
        %
        %1: Air temperature: Convert from Kelvin to Celsius
        tair(n,:,:)=tair_n-273.15;
        %
        % 2: Relative humidity: Convert from % to fraction
        rhum(n,:,:)=rhum_n/100;
        %
        % 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
        prate_n=prate_n*0.1*(24*60*60.0);
        prate_n(abs(prate_n)<1.e-4)=0;
        prate(n,:,:)=prate_n;
        %
        % 4: Net shortwave flux: [W/m^2]
        %    CROCO convention: positive downward: same as GFS
        % ?? albedo ??
        radsw_n=dradsw_n - uradsw_n;
        radsw_n(radsw_n<1.e-10)=0;
        radsw(n,:,:)=radsw_n;
        %
        % 5: Net outgoing Longwave flux:  [W/m^2]
        %    CROCO convention: positive upward (opposite to nswrs)
        %    GFS convention: positive downward --> * (-1)
        %    input: downward longwave rad. and
        %    skin temperature.
        radlw(n,:,:) = uradlw_n - dradlw_n;
        radlw_in(n,:,:)=dradlw_n;
        %
        % 6: Wind speed
        wspd(n,:,:)=sqrt(u_n.^2+v_n.^2);
        %
        % 7:  Wind vectors
        uwnd(n,:,:)=u_n;
        vwnd(n,:,:)=v_n;
        %
        % 7: Wind stress
        tx(n,:,:)=tx_n;
        ty(n,:,:)=ty_n;
        
        clear nco;
        
    else
        disp(['warning: ',gfs_grb_name,' does not exist']);
        disp('using data from the nearest available time-step');
        tair(n,:,:)=tair(n-1,:,:); % n-1 will work even in the event of multiple missing time-steps as data from n-2 will get carried over to n-1 etc. 
        rhum(n,:,:)=rhum(n-1,:,:);
        prate(n,:,:)=prate(n-1,:,:);
        radsw(n,:,:)=radsw(n-1,:,:);
        radlw_in(n,:,:)=radlw_in(n-1,:,:);
        wspd(n,:,:)=wspd(n-1,:,:);
        uwnd(n,:,:)=uwnd(n-1,:,:);
        vwnd(n,:,:)=vwnd(n-1,:,:);
        tx(n,:,:)=tx(n-1,:,:);
        ty(n,:,:)=ty(n-1,:,:);
    end
end

% check the time
check_time=abs(gfstime-gfstime_check);
if max(check_time(:))>1e-3
    disp('error in processing GFS times - debug please!')
    return
end

%
% Put the time in Yorig time
%
gfstime=gfstime-datenum(Yorig,1,1);
%
% Create the GFS output file and write everything down
%
mask(isnan(mask))=0;
write_GFS(gfs_name,Yorig,lon,lat,mask,gfstime,tx,ty,tair,rhum,prate,wspd,uwnd,vwnd,radlw,radlw_in,radsw)
%
disp('reformatting of GFS: done')
%
return




