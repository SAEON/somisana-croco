function extract_CFSR_hourly(CFSR_dir,url,vname,Y,M,lon,lat,...
                             i1min,i1max,i2min,i2max,i3min,i3max,...
                             jmin,jmax,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% Extract hourly CFSR data downloaded to your computer
% Based on extract_CFSR.m
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
%var=nan*zeros(N,Mm,Lm); % initialise variable- (time,lat,lon)
%time=nan*zeros(N,1); % initialise time datenum
%

for Dstart=1:5:26 
    
    if Dstart==26
        Dend=nmax;
    else
        Dend=Dstart+4;
    end
    
    ref_time=datenum(Y,M,Dstart);
    
    for fcst=1:6
        
        fname=[url,'/flxf',sprintf('%02d',fcst),'.gdas.',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',Dstart),'-',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',Dend),'.grb2.nc'];     
        
        var_i=read_CFSR_hourly(fname,vname,':',i1min,i1max,i2min,i2max,i3min,i3max,jmin,jmax);
        time_i=(ncread(fname,'time0'))/24.+ref_time-datenum(Yorig,1,1);
        
        % precipitation is only available as the average from the analysis
        % to the forecast hour, so we need to convert to hourly averages
        if strcmp(vname,'PRATE_L1_Avg_1')
            if fcst==1
                var_i_ave_prev=var_i;
            else
                try
                    % see FAQ "How can the individual one-hour averages be computed
                    % (https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html)
                    var_i_ave=var_i;
                    var_i=var_i_ave*(fcst)-var_i_ave_prev*(fcst-1);
                    var_i_ave_prev=var_i_ave;
                catch 
                    % we need this catch for one month 2007-12 where we don't
                    % have data written to the file for the last forecast
                    % time-step for some reason. So just repeat the last
                    % time-step as a hack around
                    var_i_ave=var_i;
                    var_i_ave=cat(1,var_i_ave,var_i_ave(end,:,:)); % here's the hack
                    var_i=var_i_ave*(fcst)-var_i_ave_prev*(fcst-1);
                    var_i_ave_prev=var_i_ave;
                end
            
            end
        end
        
        if Dstart==1 && fcst==1 % initialise output variable if first file in month
            time_all=time_i;
            var_all=var_i;
        else
            time_all=cat(1,time_all,time_i);
            var_all=cat(1,var_all,var_i);
        end
    end
end

% get the time steps in the correct order
[time_all_sort,indx]=sort(time_all);
var_all_sort=var_all(indx,:,:);
%
% Write it in a file
%
write_NCEP([CFSR_dir,vname,'_Y',num2str(Y),'M',num2str(M),'.nc'],...
	   vname,lon,lat,time_all_sort,var_all_sort,Yorig)
%
