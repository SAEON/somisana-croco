function download_GLORYS_Mydata(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
                       OGCM_dir,OGCM_prefix,path,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract a subgrid from GLORYS to get a ROMS forcing
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
    
    % Get the daily timesteps for this month
    Nt = eomday(Y,M);
    
    % first file in this month (for getting grid info)
    glorys_day=datenum(Y,M,1)-datenum(1990,1,1);
    glorysfile=[path,'algoa_',num2str(glorys_day,'%06d'),'.nc'];
    
    %
    % Find a subset of the GLORYS grid
    %
    [i1min,i1max,i2min,i2max,i3min,i3max,jrange,krange,lonT,latT,lonU,latU,lonV,latV,depth]=...
     get_GLORYS_subgrid_Mydata(glorysfile,lonmin,lonmax,latmin,latmax);
    %
    
    time=datenum(Y,M,1):datenum(Y,M,Nt);
    GLORYS_time=time-datenum(Yorig,1,1);
    trange=1:Nt;
%
% Extract GLORYS data
%
    extract_GLORYS_Mydata(OGCM_dir,OGCM_prefix,path,Y,M,...
                 lonT,latT,lonU,latU,lonV,latV,depth,GLORYS_time,...
                 trange,krange,jrange,...
                 i1min,i1max,i2min,i2max,i3min,i3max,...
                 Yorig)
  end
end
return
