function download_CFSR_hourly(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
    NCEP_dir,Yorig,My_NCEP_dir,NCEP_version)
%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% script is based on download_CFSR.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data from global netcdf format saved on my computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Illig, 2010, from download_NCEP
%  Updated    January 2016 (E. Cadier and S. Illig)
%  Updated    May 2017 (S. Illig)
%

%
% Definitions of names and directories
%
disp(['===================='])
disp(['Direct CONVERT Procedure'])
disp(['===================='])

disp(['CFSR directory:'])
ncep_url=My_NCEP_dir

%
% Definitions of names and directories
%
if NCEP_version==4
%
  vnames={'TMP_L1' ...      % surface temp. [k]
      'DLWRF_L1' ...   % surface downward long wave flux [w/m^2]
      'ULWRF_L1' ...   % surface upward long wave flux [w/m^2]
      'DSWRF_L1' ...   % surface downward solar radiation flux [w/m^2]
      'USWRF_L1' ...   % surface upward solar radiation flux [w/m^2]
      'PRATE_L1_Avg_1' ...   % surface precipitation rate, average from analysis to forecast tstep [kg/m^2/s]
      'U_GRD_L103' ...    % 10 m u wind [m/s]
      'V_GRD_L103' ...    % 10 m v wind [m/s]
      'TMP_L103' ...      % temperature at 2m above surface [k]
      'SPF_H_L103'};       % 2 m specific humidity [kg/kg]
else
  error('Wrong NCEP version')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common OpenDAP FTP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' '])
disp(['Get CFSR data from ',num2str(Ymin),' to ',num2str(Ymax)])
disp(['From ',ncep_url]);
disp([' '])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Making output data directory ',NCEP_dir])
eval(['!mkdir ',NCEP_dir])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find a subset of the CFSR grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifile=[ncep_url,'/flxf01.gdas.',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'01-',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'05.grb2.nc'];
[i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax,lon,lat]=...
 get_CFSR_hourly_subgrid(ifile,lonmin,lonmax,latmin,latmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global loop on variable names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(vnames)
    disp(['=========================='])
    disp(['VNAME IS ',char(vnames(k))]);
    disp(['=========================='])
    vname=char(vnames(k));
    if k==1 % Dealing with the mask (% WHY IS THE LAND MASK NEEDED FOR ATMOSPHERIC FORCING?)
        %
        disp(['==========================']);
        disp(['Get the Land Mask tindex = 1']);
        disp(['Get the Land Mask by using SSTK']);
        disp([' '])
        %
        ifile=[ncep_url,'/flxf01.gdas.',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'01-',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'05.grb2.nc'];
        
        var=read_CFSR_hourly(ifile,vname,1,i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
        var(var == min(min(var)))=NaN;
        mask=var-var;
        % just use the ECMWF method for writing a mask file
        write_ECMWF_Mask([NCEP_dir,'land_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'],...
            'land',lon,lat,mask)
        disp([' '])
        clear var;
        %
    end
  
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

       extract_CFSR_hourly(NCEP_dir,ncep_url,...
                    char(vnames(k)),Y,M,lon,lat,...
                    i1min,i1max,i2min,i2max,i3min,i3max,...
                    jmin,jmax,Yorig)
      end % end loop month
    end % end loop year
end % loop k

return
