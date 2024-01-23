function download_GLORYS_Mydata(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
                       OGCM_dir,OGCM_prefix,path,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract a subgrid from GLORYS files to get a ROMS forcing
% Store that into monthly files.
% Take care of the Greenwitch Meridian.
% 
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
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
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' '])
disp(['Get data from Y',num2str(Ymin),'M',num2str(Mmin),...
      ' to Y',num2str(Ymax),'M',num2str(Mmax)])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])
%

%
% Create the directory
%
disp(['Making output data directory ',OGCM_dir])
eval(['!mkdir ',OGCM_dir])
%
% Start 
%
disp(['Process the dataset: ',path])

%
% Loop on the years
%
for Y=Ymin:Ymax
  disp(['Processing year: ',num2str(Y)])
%
% Loop on the months

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
    
    glorysfile=[path,num2str(Y),'_',num2str(M,'%02.f'),'.nc'];
    
    %
    % Find a subset of the HYCOM grid
    %
    [i1min,i1max,i2min,i2max,i3min,i3max,jrange,krange,lon,lat,depth]=...
     get_GLORYS_subgrid_Mydata(glorysfile,lonmin,lonmax,latmin,latmax);
    %
    % Get GLORYS time 
    %
    nc=netcdf(glorysfile);
    time=nc{'time'}(:);
    nc(close);
    time=time/24; % hours to days
    Nt=length(time);
    time_origin=ncreadatt(glorysfile,'time','units');
    time_origin=strsplit(time_origin);
    time_origin=char(strcat(time_origin(3)," ",time_origin(4)));
    time_origin=datenum(time_origin);
    % Transform time to days since hycom time origin
    time=time+time_origin;
    % Transform it into Yorig time (i.e days since Yorig-01-01)
    GLORYS_time=time-datenum(Yorig,1,1);
    
    trange=1:Nt;
%
% Extract GLORYS data
%
    extract_GLORYS_Mydata(OGCM_dir,OGCM_prefix,glorysfile,Y,M,...
                 lon,lat,depth,GLORYS_time,...
                 trange,krange,jrange,...
                 i1min,i1max,i2min,i2max,i3min,i3max,...
                 Yorig)
  end
end
return
