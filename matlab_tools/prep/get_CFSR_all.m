function [var_all_sort, time_all_sort] =get_CFSR_all(url,vname,datetime_start,datetime_end)

% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% this function returns already downloaded hourly CFSR data over the
% specified time period (no flipping of lat dimension)

year_start=year(datetime_start);
year_end=year(datetime_end);
month_start=month(datetime_start);
month_end=month(datetime_end);

for Y=year_start:year_end
    
    if Y==year_start
        M_start=month_start;
    else
        M_start=1;
    end
    if Y==year_end
        M_end=month_end;
    else
        M_end=12;
    end
    
    for M=M_start:M_end
        %
        nmax=daysinmonth(Y,M);
        %
        
        for Dstart=1:5:26
            
            if Dstart==26
                Dend=nmax;
            else
                Dend=Dstart+4;
            end
            
            ref_time=datenum(Y,M,Dstart);
            
            for fcst=1:6
                
                fname=[url,'flxf',sprintf('%02d',fcst),'.gdas.',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',Dstart),'-',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',Dend),'.grb2.nc'];
                
                nc=netcdf(fname);
                var_i=nc{vname}(:);
                time_i=(ncread(fname,'time0'))/24.+ref_time; %need to change to time0 if data from different averaging periods in same file
                close(nc);
                
                % precipitation is only available as the average from the analysis
                % to the forecast hour, so we need to convert to hourly averages
                if strcmp(vname,'PRATE_L1_Avg_1')
                    if fcst==1
                        var_i_ave_prev=var_i;
                    else
                        % see FAQ "How can the individual one-hour averages be computed
                        % (https://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html)
                        var_i_ave=var_i;
                        var_i=var_i_ave*(fcst)-var_i_ave_prev*(fcst-1);
                        var_i_ave_prev=var_i_ave;
                    end
                end
                
                if Y==year_start && M==month_start && Dstart==1 && fcst==1 % initialise output variable if first file
                    time_all=time_i;
                    var_all=var_i;
                else
                    time_all=cat(1,time_all,time_i);
                    var_all=cat(1,var_all,var_i);
                end
                
            end
        end
        
    end
end

% get the time steps in the correct order
[time_all_sort,indx]=sort(time_all);
var_all_sort=var_all(indx,:,:);

% only keep the timesteps between specified dates
datetime_all_sort=datetime(time_all_sort, 'ConvertFrom', 'datenum');
trange=find(isbetween(datetime_all_sort,datetime_start,datetime_end));
datetime_all_sort=datetime_all_sort(trange);
time_all_sort=time_all_sort(trange);
var_all_sort=var_all_sort(trange,:,:);

%
return