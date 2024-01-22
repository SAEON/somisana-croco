function [time,data] = get_ts_indx(dirname,prefix,I_indx,J_indx,K_indx,vname,time_origin,time_start,time_end)
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts
% extracts at defined I, J and K indices specified by the user

% get the grid info from the first file (common to all)
year_now=year(datetime(time_start,'ConvertFrom','datenum'));
month_now=month(datetime(time_start,'ConvertFrom','datenum'));

% extract the data by looping through the monthly files. This gets all data
% in the monthly files. We cut off the uneeded data afterward
date_eof = time_start; % initialise the time at the end of the current file
count=1;
while date_eof < time_end
    
    crocofile=[dirname,prefix,'_Y',num2str(year_now),'M',num2str(month_now),'.nc'];
    nc=netcdf(crocofile);
    
    % time
    time_croco=nc{'time'}(:);
    time_croco=time_croco/3600/24; % days
    time_croco=time_croco+time_origin;

    % read the variable
    var=squeeze(nc{vname}(:,K_indx,J_indx,I_indx));

    if count==1
        time=time_croco;
        data=var;
    else
        time=cat(1,time,time_croco);
        data=cat(1,data,var);
    end
    
    count=count+1;
    month_now=month_now+1;
    if month_now==13
        year_now=year_now+1;
        month_now=1;
    end
        
    % check if we need to keep going
    date_eof=time(end);
    
    close(nc);
end

% now subset the data using the start and end times
indx=find(time>=time_start&time<=time_end);
time=time(indx);
data=data(indx);

end