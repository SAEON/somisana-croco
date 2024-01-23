function datenum_rounded = datenum_round_off(datenum_in, time_unit_string)
% Rounds off a datenum to the nearest second, minute, hour, whatever.
%
% datenum_rounded = datenum_round_off(datenum_in, time_unit_string)
%
% Inputs:
%   datenum_in          Date/time in Matlab datenum format
%   time_unit_string    String like 'minute', 'second', etc.
%
% Output:
%   datenum_rounded     Input value rounded off

% Kevin J. Delaney
% January 13, 2010

datenum_rounded = [];

if ~exist('datenum_in', 'var')
    help(mfilename);
    return
end

if isempty(datenum_in) || ~isnumeric(datenum_in)
    errordlg('Input "datenum_in" is empty or non-numeric.', mfilename);
    return
end

if ~exist('time_unit_string', 'var') || ...
   isempty(time_unit_string) || ...
   ~ischar(time_unit_string)
    errordlg('Input "time_unit_string" is missing, empty or non-char.', ...
        mfilename);
    return
end

switch lower(time_unit_string)
    case {'second', 'seconds', 'sec', 's'}
        time_unit = units(1, 'second', 'day');
        
    case {'minute', 'minutes', 'min', 'm'}
        time_unit = units(1, 'minute', 'day');

    case {'hour', 'hours', 'hr', 'h'}
        time_unit = units(1, 'hour', 'day');

    case {'day', 'days', 'd'}
        time_unit = 1;

    case {'month', 'months', 'mon'}
        %   Break up the input into year, month, day components.
        [yr, mon, day, hr, min, sec] = datevec(datenum_in);
        
        %   Are we halfway through the month?
        datenum_start_of_this_month = datenum(yr, mon, ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        datenum_start_of_next_month = datenum(yr, mon + 1, ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        round_up_to_next_month_syndrome = (datenum_start_of_next_month - datenum_in) < (datenum_in - datenum_start_of_this_month);
        boost_vector = zeros(size(datenum_in));
        boost_vector(round_up_to_next_month_syndrome) = 1;
        datenum_rounded = datenum(yr, mon + boost_vector, ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        return

    case {'year', 'years', 'yr'}
        %   Break up the input into year, month, day components.
        [yr, mon, day, hr, min, sec] = datevec(datenum_in);
        
        %   Are we halfway through the year?
        datenum_start_of_this_year = datenum(yr, ones(size(mon)), ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        datenum_start_of_next_year = datenum(yr + 1, ones(size(mon)), ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        round_up_to_next_year_syndrome = (datenum_start_of_next_year - datenum_in) < (datenum_in - datenum_start_of_this_year);
        boost_vector = zeros(size(datenum_in));
        boost_vector(round_up_to_next_year_syndrome) = 1;
        datenum_rounded = datenum(yr + boost_vector, ones(size(mon)), ones(size(day)), zeros(size(hr)), zeros(size(min)), zeros(size(sec)));
        return
        
    otherwise
        errordlg(['Unknown time unit: "', time_unit_string, '".'], mfilename);
        return
end

datenum_rounded = time_unit * round(datenum_in / time_unit);