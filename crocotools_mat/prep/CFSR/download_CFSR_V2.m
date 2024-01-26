function download_CFSR_V2(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
                       NCEP_dir,NCEP_version,Yorig,Get_My_Data,My_NCEP_dir)
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% This is a quick job to use CFSR V2 data already processed by Gildas
% Cambon... this script just interpolates the data to CROCO input files
%
%clear all
%close all
if nargin < 1
  Ymin=2006;
  Ymax=2006;
  Yorig=1900;
  Mmin=6;
  Mmax=6;
  lonmin=16;
  lonmax=19;
  latmin=-34;
  latmax=-32;
  NCEP_dir='DATA/CFSR/';
  NCEP_version=3;
  Get_My_Data=0;
  My_NCEP_dir='../CFSR/';
end
%
% Definitions of names and directories
%
if NCEP_version==3
%
% Definitions of names and directories for NCEP1
%

  ncep_url='/home/gilesfearon/phd/data/CFSR/6_hourly_V2/';
  vnames={'Temperature_height_above_ground' ...      % 2 m temp. [k]
      'Downward_Long-Wave_Rad_Flux' ...   % surface downward long wave flux [w/m^2]
      'Upward_Long-Wave_Rad_Flux_surface'   ...   % surface upward long wave flux [w/m^2] 'Temperature_surface' ...     % surface temp. [k]
      'Downward_Short-Wave_Rad_Flux_surface' ...   % surface downward solar radiation flux [w/m^2]
      'Upward_Short-Wave_Rad_Flux_surface' ...   % surface upward solar radiation flux [w/m^2]
      'Precipitation_rate' ...   % surface precipitation rate [kg/m^2/s]
      'U-component_of_wind' ...    % 10 m u wind [m/s]
      'V-component_of_wind' ...    % 10 m v wind [m/s]
      'Specific_humidity'};       % 2 m specific humidity [kg/kg]

else
  error('Wrong NCEP version')
end

disp([' '])
disp(['Get CFSR data from ',num2str(Ymin),' to ',num2str(Ymax)])
disp(['From ',ncep_url]);
disp([' '])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])
%
% Create the directory
%
disp(['Making output data directory ',NCEP_dir])
eval(['!mkdir ',NCEP_dir])
% End Common OpenDap FTP 
%-----------------------

%
% Global loop on variable names
%
for k=1:length(vnames)
  disp(['=========================='])
  disp(['VNAME IS ',char(vnames(k))]);
  disp(['=========================='])

  if k==1
    disp(['=========================='])
    disp(['GET  SUBGRID time only k=1']);
    disp(['USE VARIABLE: ',char(vnames(k))  ])
    disp(['=========================='])
    
    fname=[ncep_url,char(vnames(k)),'_Y',num2str(Ymin),'.nc'];
%
    [imin,imax,jmin,jmax,lon,lat]=...
    get_NCEP_grid_V2(fname,...
    lonmin,lonmax,latmin,latmax);

  end %k==1
%
% Loop on the years
%
    for Y=Ymin:Ymax
      disp(['=========================='])
      disp(['Processing year: ',num2str(Y)])
      disp(['=========================='])
      
      fname=[ncep_url,char(vnames(k)),'_Y',num2str(Y),'.nc'];
      time=(ncread(fname,'time'))+datenum(1900,1,1);
      %time_datestr=datestr(time);
     
      nc=netcdf(fname);
      
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
          ndays = eomday(Y,M);
          start_datenum=datenum(Y,M,1);
          end_datenum=datenum(Y,M,ndays,18,0,0);
          
          time_sub_indx=find(time>=start_datenum & time<=end_datenum);
          time_sub=time(time_sub_indx);
          time_sub=time_sub-datenum(Yorig,1,1);
          var_sub=nc{char(vnames(k))}(time_sub_indx,jmin:jmax,imin:imax);
          
          write_NCEP([NCEP_dir,char(vnames(k)),'_Y',num2str(Y),'M',num2str(M),'.nc'],...
	      char(vnames(k)),lon,lat,time_sub,var_sub,Yorig)
        
      end % end loop month
      close(nc)
    end % end loop year
end % loop k
%
return
