function extract_CFSRV2_hourly(CFSR_dir,url,vname,Y,M,lon,lat,...
                             i1min,i1max,i2min,i2max,i3min,i3max,...
                             jmin,jmax,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract hourly CFSR data downloaded to your computer
% Based on extract_CFSR.m as part of crocotools.m
%
% Write it in a local file
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
%  Copyright (c) 2011 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Get the variable name
%
disp(['Get ',vname,' for year ',num2str(Y),...
      ' - month ',num2str(M)])
%
% Get the number of days in the month
%      
nmax=daysinmonth(Y,M);
Lm=length(lon);
Mm=length(lat);
N=24*nmax;
%
for D=1:nmax
    
    ref_time=datenum(Y,M,D);
    fname=[url,'cdas1.',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',D),'.sfluxgrbf.grb2.nc'];
    
    nc=netcdf(fname);
    
    var_i=read_CFSR_hourly(fname,vname,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
    time_i=(ncread(fname,'time0'))/24.+ref_time-datenum(Yorig,1,1);
    
    forecast_hour=nc{'forecast_hour0'}(:);
    % subset time and var to remove unwanted forecast hours if they
    % exist in the file
    indx=find(forecast_hour<=6);
    time_i=time_i(indx);
    var_i=var_i(indx,:,:);
    
    if D==1 % initialise output variable if first file
        time_all=time_i;
        var_all=var_i;
    else
        time_all=cat(1,time_all,time_i);
        var_all=cat(1,var_all,var_i);
    end
    close(nc)
end

%
% Write it in a file
%
write_NCEP([CFSR_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc'],...
	   vname,lon,lat,time_all,var_all,Yorig)
%
