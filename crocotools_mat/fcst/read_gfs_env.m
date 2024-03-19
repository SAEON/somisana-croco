function [run_date,delta_days,hdays,fdays]=read_gfs_env(gfs_env)

fileID = fopen(gfs_env, 'r');
% Read the file line by line
while true 
    thisline = fgetl(fileID); 
    % Check if we've reached the end of the file
    if ~ischar(thisline)
        break;
    end
    % Split the line into variable name and value
    parts = strsplit(thisline, '=');
    if length(parts) == 2
        varName = strtrim(parts{1});
        varValue = strtrim(parts{2});
        % **3. Process each variable**
        switch varName
            case 'RUN_DATE'
                % Convert datestring to datenum
                run_date = datenum(varValue, 'yyyy-mm-dd HH');  
            case 'DELTA_DAYS_GFS'
                delta_days = str2double(varValue); 
            case 'HDAYS'
                hdays = str2double(varValue);
            case 'FDAYS'
                fdays = str2double(varValue);
        end
    end
end
fclose(fileID);
