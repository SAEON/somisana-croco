function download_ECMWF_3_hourly(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
    ECMWF_dir,Yorig,My_ECMWF_dir)
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract data from global netcdf format saved on my computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Illig, 2010, from download_NCEP
%  Updated    January 2016 (E. Cadier and S. Illig)
%  Updated    May 2017 (S. Illig)
%
%  Edited from download_ECMWF as part of roms-tools to include 3 hourly forecast
%

if nargin < 1
    Ymin=2007;
    Ymax=2008;
    Yorig=1900;
    Mmin=12;
    Mmax=1;
    lonmin=10;
    lonmax=30;
    latmin=-40;
    latmax=-20;
    ECMWF_dir='DATA/NCEP_Peru/';
    My_ECMWF_dir='../NCEP_REA1/';
end
%
% Definitions of names and directories
%
disp(['===================='])
disp(['Direct CONVERT Procedure'])
disp(['===================='])

disp(['ECMWF directory:'])
ecmwf_url=My_ECMWF_dir
catalog={'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...    
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...
    'EI_ecmwf_' ...    
    ''};
vnames={'SSTK'...  % surface land-sea mask [1=land; 0=sea]
    'T2M' ...      % 2 m temp. [k]
    'SSTK' ...     % surface temp. [k]
    'U10M' ...     % 10 m u wind [m/s]
    'V10M' ...     % 10 m v wind [m/s]
    'Q' ...        % 2 m specific humidity [kg/kg]
    'STR' ...      % surface thermal radiation [w/m^2.s] -> [w/m^2]
    'STRD' ...     % surface downward thermal radiation [w/m^2.s] -> [w/m^2]
    'SSR' ...      % surface solar radiation [w/m^2.s] -> [w/m^2]
    'TP' ...       % surface precipitation rate [m]->[kg/m^2/s]
    'EWSS' ...     % east-west surface stress [N m-2 s] -> [N m-2]
    'NSSS' ...     % north-south surface stress [N m-2 s] -> [N m-2]
    };          % 2 m specific humidity [kg/kg]
fnames={'sst'...  % surface land-sea mask [1=land; 0=sea]
    't2m' ...      % 2 m temp. [k]
    'sst' ...     % surface temp. [k]
    'u10' ...     % 10 m u wind [m/s]
    'v10' ...     % 10 m v wind [m/s]
    'q' ...        % 2 m specific humidity [kg/kg]
    'str' ...      % surface thermal radiation [w/m^2.s] -> [w/m^2]
    'strd' ...     % surface downward thermal radiation [w/m^2.s] -> [w/m^2]
    'ssr' ...      % surface solar radiation [w/m^2.s] -> [w/m^2]
    'tp' ...       % surface precipitation rate [m]->[kg/m^2/s]
    'ewss' ...     % east-west surface stress [N m-2 s] -> [N m-2]
    'nsss' ...     % north-south surface stress [N m-2 s] -> [N m-2]
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common OpenDAP FTP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' '])
disp(['Get ECMWF data from ',num2str(Ymin),' to ',num2str(Ymax)])
disp(['From ',ecmwf_url]);
disp([' '])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Making output data directory ',ECMWF_dir])
eval(['!mkdir ',ECMWF_dir])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find a subset of the ECMWF grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ifile=[char(ecmwf_url),char(catalog(1)),char(vnames(1)),'_',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'.nc'];
[i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax,lon,lat]=...
 get_ECMWF_subgrid(ifile,lonmin,lonmax,latmin,latmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global loop on variable names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:length(vnames)
    disp(['=========================='])
    disp(['VNAME IS ',char(vnames(k))]);
    disp(['=========================='])
    %
    % 
    vname=char(vnames(k));
    fname=char(fnames(k));
    if k==1 % Dealing with the mask (using SST)
        %
        disp(['==========================']);
        disp(['Get the Land Mask tindex = 1']);
        disp(['Get the Land Mask by using SSTK']);
        disp([' '])
        %
        ifile=[char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Ymin),sprintf('%02d',Mmin),'.nc'];
        var=extract_ECMWF(ifile,fname,Ymin,Mmin,1,i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
        var(var == min(min(var)))=NaN;
        mask=var-var;
        write_ECMWF_Mask([ECMWF_dir,'land_Y',num2str(Ymin),'M',num2str(Mmin),'.nc'],...
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %
        for M=mo_min:mo_max
            %
            if k<=6   % t2m, sst, u10, v10, q
                
                file =[ECMWF_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc'];
                ifile=[char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'.nc'];
                ifile3=[char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step3.nc'];
                %
                var0=extract_ECMWF(ifile,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
                var3=extract_ECMWF(ifile3,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
                time0=(ncread(ifile,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                time3=(ncread(ifile3,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                %
                % subset to timesteps for this month only
                time0=time0(5:end-4);
                var0=var0(5:end-4,:,:);
                time3=time3(5:end-4);
                var3=var3(5:end-4,:,:);
                %
                % slot the forecast timesteps in between the 6 hourly 
                % analysis timesteps
                var(1:2:(length(time0))*2-1,:,:)=var0; clear var0;
                var(2:2:(length(time3))*2,:,:)=var3; clear var3;          
                time(1:2:(length(time0))*2-1)=time0; clear time0;
                time(2:2:(length(time3))*2)=time3; clear time3;
                
                write_ECMWF(file,vname,lon,lat,time,var,Yorig);
                clear time;clear var;
                clear file; clear ifile;
                disp([' ']);
            %
            else  % str, strd, ssr, tp, ewss, nsss
        
                file=[ECMWF_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc'];
                %
                ifile3  = [char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step3.nc'];
                ifile6  = [char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step6.nc'];
                ifile9  = [char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step9.nc'];
                ifile12  = [char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step12.nc'];
                ifile15 = [char(ecmwf_url),char(catalog(k)),vname,'_',sprintf('%04d',Y),sprintf('%02d',M),'_step15.nc'];
                %
                div   = (86400./8.); % 8 three hour timesteps in a day
                if strcmp(vname,'TP') 
                    div=div/1000.;
                end
                if strcmp(vname,'STR') 
                    div=-1.*div;
                end 
                %
                var3  = extract_ECMWF(ifile3,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax)/div;  %divide by 3h;
                var6  = extract_ECMWF(ifile6,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax)/div;  %divide by 3h;
                var9  = extract_ECMWF(ifile9,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax)/div;  %divide by 3h;
                var12  = extract_ECMWF(ifile12,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax)/div;  %divide by 3h;
                %var15 = extract_ECMWF(ifile15,fname,Y,M,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax)/div; %divide by 3h;
                time3=(ncread(ifile3,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                time6=(ncread(ifile6,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                time9=(ncread(ifile9,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                time12=(ncread(ifile12,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                %time15=(ncread(ifile15,'time'))/24.+datenum(1900,1,1)-datenum(Yorig,1,1);
                %
                
                var(1:4:length(time3)*4-3,:,:)=var3;
                var(2:4:length(time3)*4-2,:,:)=var6-var3; clear var3;
                var(3:4:length(time3)*4-1,:,:)=var9-var6; clear var6;
                var(4:4:length(time3)*4,:,:)=var12-var9; clear var9; clear var12;
                time(1:4:length(time3)*4-3)=time3-1/16.; % shift time back by 1.5 hours (1/16 days) for the rate over the first 3 hour accumulation
                time(2:4:length(time3)*4-2)=(time6+time3)/2.; 
                time(3:4:length(time3)*4-1)=(time9+time6)/2.; clear time6;
                time(4:4:length(time3)*4)=(time12+time9)/2.; clear time3; clear time9; clear time12;

                % using 3 hourly accumulated data as above means that we
                % get rates every 3 hours, but starting from the 1.5 hr
                % timestep. I prefer to do it this way rather than using accumulated
                % data over a 6 hour period, which gets you the correct
                % timestep, but the rate is averaged over a longer period
                % So we need to interpolate to get back to the 3 hour timestep
                % to be consistent with variables 1-6
                var=(var(1:end-1,:,:)+var(2:end,:,:))/2.;
                time=(time(1:end-1)+time(2:end))/2.;
                
                time2=time(8:end-8);clear time;
                var2=var(8:end-8,:,:);clear var;
                %
                write_ECMWF(file,vname,lon,lat,time2,var2,Yorig);
                clear time2; clear var2;clear div;
                clear file; clear ifile3; clear ifile9; clear ifile15;
                disp([' ']);
                %
            end % end if
        end % end loop month
    end % end loop year
end % loop k
%
return
