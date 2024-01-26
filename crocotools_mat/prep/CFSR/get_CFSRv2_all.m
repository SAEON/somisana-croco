function [var_all, time_all] =get_CFSRv2_all(url,vname,datetime_start,datetime_end)

% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% this function returns already downloaded hourly CFSR V2 data over the
% specified time period (no flipping of lat dimension)

year_start=year(datetime_start);
year_end=year(datetime_end);
month_start=month(datetime_start);
month_end=month(datetime_end);
day_start=day(datetime_start);
day_end=day(datetime_end);

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
        if M==month_start
        D_start=day_start;
        else
            D_start=1;
        end
        if M==month_end
            D_end=day_end;
        else
            D_end=nmax;
        end
        
        for D=D_start:D_end
            
            ref_time=datenum(Y,M,D);
            fname=[url,'cdas1.',sprintf('%04d',Y),sprintf('%02d',M),sprintf('%02d',D),'.sfluxgrbf.grb2.nc'];
                
            nc=netcdf(fname);
            var_i=nc{vname}(:);
            time_i=(ncread(fname,'time0'))/24.+ref_time;
            forecast_hour=nc{'forecast_hour0'}(:);
            close(nc);
            % subset time and var to remove unwanted forecast hours if they
            % exist in the file
            indx=find(forecast_hour<=6);
            time_i=time_i(indx);
            var_i=var_i(indx,:,:);
            
            if Y==year_start && M==month_start && D==day_start % initialise output variable if first file
                time_all=time_i;
                var_all=var_i;
            else
                time_all=cat(1,time_all,time_i);
                var_all=cat(1,var_all,var_i);
            end
        end
    end
end

% only keep the timesteps between specified dates
% (only relevent if start and end time HH:MM:SS are not 00:00:00)
datetime_all=datetime(time_all, 'ConvertFrom', 'datenum');
trange=find(isbetween(datetime_all,datetime_start,datetime_end));
datetime_all=datetime_all(trange);
time_all=time_all(trange);
var_all=var_all(trange,:,:);

%
return