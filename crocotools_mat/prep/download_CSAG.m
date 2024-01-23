function download_CSAG(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
    CSAG_dir,Yorig,My_CSAG_dir)
%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% script is based on download_CFSR.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data from subset of netcdf's found on terra
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%
% Definitions of names and directories
%
disp(['===================='])
disp(['Direct CONVERT Procedure'])
disp(['===================='])

disp(['CSAG directory:'])
csag_url=My_CSAG_dir

%
% Definitions of names and directories
%

vnames={'rlds' ...   % downward long wave flux at ground surface [w.m^-2]
    'rsds' ...   % downward short wave flux at ground surface [w.m^-2]
    'rain' ...   % rainfall [mm]
    'u10' ...    % 10 m u wind [m/s]
    'v10' ...    % 10 m v wind [m/s]
    'tas' ...      % temperature at 2m above surface [k]
    'huss'};       % 2 m specific humidity [kg/kg]

% NOTE: WE DON"T HAVE UPWARD FLUXES FROM CSAG- NEED TO ESTIMATE THESE IN
% INTERP_CSAG.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common OpenDAP FTP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' '])
disp(['Get CFSR data from ',num2str(Ymin),' to ',num2str(Ymax)])
disp(['From ',csag_url]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Making output data directory ',CSAG_dir])
eval(['!mkdir ',CSAG_dir])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the CSAG grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
nc=netcdf([csag_url,'/wrfout_d03_wind_10m_regridded_sub.nc']); % any file will do - grid same for all
lon=squeeze(nc{'lon'}(:));
lat=squeeze(nc{'lat'}(:));
close(nc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global loop on variable names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(vnames)
    disp(['=========================='])
    disp(['VNAME IS ',char(vnames(k))]);
    disp(['=========================='])
    vname=char(vnames(k));
    if vname=='u10' || vname=='v10'
        fname=[csag_url,'/wrfout_d03_wind_10m_regridded_sub.nc'];
    else
        fname=[csag_url,'/wrfout_d03_',vname,'wind_10m_regridded_sub.nc'];
    end
    %
    disp('==========================')
    disp(['getting time array for: ',vname])
    disp('==========================')
    time=ncread(fname,'time');
    time=time/24; % days since reference time
    ref_time=datenum(1989,1,1);
    time=time+ref_time; % time as matlab datenum
    %
    % noticed some duplicate time-steps so remove these
    [time, time_indx, ~] = unique(time);
    datestrings=datestr(time);
    datetimes=datetime(datestrings);
    time=time-datenum(Yorig,1,1); % time is now in days since 1-Jan-Yorig
    
%    if k==1 % Dealing with the mask (% WHY IS THE LAND MASK NEEDED FOR ATMOSPHERIC FORCING?)
%        %
%         disp(['==========================']);
%         disp(['Get the Land Mask tindex = 1']);
%         disp(['Get the Land Mask by using SSTK']);
%         disp([' '])
%         %
%         ifile=[csag_url,'/flxf01.gdas.',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'01-',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'05.grb2.nc'];
%         
%         var=read_CFSR_hourly(ifile,vname,1,i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
%         var(var == min(min(var)))=NaN;
%         mask=var-var;
%         % just use the ECMWF method for writing a mask file
%         write_ECMWF_Mask([CSAG_dir,'land_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'],...
%             'land',lon,lat,mask)
%         disp([' '])
%         clear var;
%        %
%    end
    
    % get time for this variable
  
%
% Loop on the years
%
    for Y=Ymin:Ymax
      disp(['=========================='])
      disp(['Processing year: ',num2str(Y)])
      disp(['=========================='])
%
% Loop on the months
%
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
        disp(['  Processing month: ',num2str(M)])
        %
        % Get the time subset for this month
        %
        % specify start time
        datenum_start=datenum(Y,M,1,0,0,0);
        datestr_start=datestr(datenum_start,0);
        datetime_start=datetime(datestr_start);
        % specify end time
        nmax=daysinmonth(Y,M);
        datenum_end=datenum(Y,M,nmax,23,0,0);
        datestr_end=datestr(datenum_end,0);
        datetime_end=datetime(datestr_end);
        % get the time indices for this month
        var_indx=find(isbetween(datetimes,datetime_start,datetime_end));
        %
        % Get the data for the variable
        %
        nc=netcdf(fname);
        var_csag=nc{vname}(time_indx(var_indx),:,:); % using time_indx(var_indx) gets the time indices of the original data
        var_missing = ncreadatt(fname,vname,'_FillValue');
        var_csag(var_csag==var_missing)=NaN;
        close(nc)
        %
        % Write it in a file (write_NCEP.m should work fine for this)
        %
        write_NCEP([CSAG_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc'],...
            vname,lon,lat,time,var_csag,Yorig)
        
      end % end loop month
    end % end loop year
end % loop k

return
