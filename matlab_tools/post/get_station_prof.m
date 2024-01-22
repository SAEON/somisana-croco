function [time,data_prof,depth_prof]=get_station_prof(dirname,sta_no,var,origin_datenum,time_start,time_end)
%
% script created by Giles Fearon, adapted from existing CROCOTOOLS scripts

year_now=year(datetime(time_start,'ConvertFrom','datenum'));
month_now=month(datetime(time_start,'ConvertFrom','datenum'));

% extract the data by looping through the monthly files. This gets all data
% in the monthly files. We cut off the uneeded data afterward
date_eof = time_start; % initialise the time at the end of the current file
count=1;
while date_eof <= time_end
    
    fname=[dirname,'sta_Y',num2str(year_now),'M',num2str(month_now),'.nc'];
    
    % get time
    sta_time=ncread(fname,'scrum_time');
    sta_time=sta_time/3600./24.; % seconds to days
    sta_time=sta_time+origin_datenum; % matlab date num
    sta_time=sta_time(2:end);  % start on 2nd tstep as first tstep contains rubbish data
    %sta_datestr=datestr(sta_time);
    
    % get station depth
    sta_depth=ncread(fname,'depth');
    sta_depth=squeeze(sta_depth(:,sta_no,2:end));
    sta_depth=sta_depth';
    
    % get station variable
    sta_var=ncread(fname,var);
    sta_var=squeeze(sta_var(:,sta_no,2:end));
    sta_var=sta_var';
    
%     sta_lon=ncread(fname,'lon');
%     sta_lon=sta_lon(sta_no);
%     sta_lat=ncread(fname,'lat');
%     sta_lat=sta_lat(sta_no);
%     disp(['getting station coords: ',num2str(sta_lon),', ',num2str(sta_lat)]);
%     
    if count==1
        time=sta_time;
        data_prof=sta_var;
        depth_prof=sta_depth;
    else
        time=cat(1,time,sta_time);
        data_prof=cat(1,data_prof,sta_var);
        depth_prof=cat(1,depth_prof,sta_depth);
    end
    
    count=count+1;
    month_now=month_now+1;
    if month_now==13
        year_now=year_now+1;
        month_now=1;
    end
        
    % check if we need to keep going
    date_eof=sta_time(end);

end

% now subset the data using the start and end times
indx=find(time>=time_start&time<=time_end);
time=time(indx);
data_prof=data_prof(indx,:);
depth_prof=depth_prof(indx,:);

disp(['model bottom layer depth at station = ',num2str(mean(depth_prof(:,1))),' m']);

return
