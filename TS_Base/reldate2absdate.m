function [tme_out, num_out] = reldate2absdate(time_in, unit_in);
% The function transforms relative dates (as e.g. from CF-conform
% netCDF-files) into absolute date-times. 
%--------------------------------------------------------------------------
% INPUT:
% - time_in     Vector (or scalar) with relative time-values
% - unit_in     String which contains the unit of the time-values. Usually,
%               this is something like "hours since 1900-01-01 12:00:00". 
%               If there is a difference between the reference date and
%               UTC, this is generally considered by adding the offset at
%               the end of the time unit (i.e. "hours since 1900-01-01
%               12:00:00 -0400"). 
%--------------------------------------------------------------------------
% OUTPUT:
% - tme_out     Matrix with absolute dates. The matrix has six rows, which
%               contain the elements YYYY MM DD HH MM SS. 
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         December 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: dtevec.m
%--------------------------------------------------------------------------

% Read the unit string from the time variable
tme_unit = textscan(unit_in, '%s %s %s %s %s');

if ~strcmp(tme_unit{1}, 'days') & ~strcmp(tme_unit{1}, 'hours') & ...
        ~strcmp(tme_unit{1}, 'seconds')
    error('Unknown time unit')
else
    base_unit = tme_unit{1};
end

% Convert the third element into a year-month-day datetime
if ~isempty(tme_unit{3})
    ymd      = datetime(tme_unit{3}, 'InputFormat', 'yyyy-MM-dd');
else
    error('Could not read yyyy-MM-dd!')
end

% Convert the fourth element into a hours-minutes-seconds datetime
if ~isempty(tme_unit{4})
    hms      = datetime(tme_unit{4}, 'InputFormat', 'HH:mm:ss');
else
    error('Could not read HH-MM-SS!')
end

% Offset to UTC
if ~isempty(tme_unit{5})
    tmp     = length(tme_unit{5}{:});
    utc_off = str2num(tme_unit{5}{:});
    
    if tmp == 5
        utc_off = utc_off/100;
    end
    hms     = hms + hours(utc_off);
end

% Combine YMD with HMS
ref_dte     = ymd + timeofday(hms);
% Compute the numeric date of this reference date
ref_dte_num = double(datenum(ref_dte));
% As relative time-vectors are usually provided in units of "seconds 
% since", "hours since", or "days since", we can simply add the original 
% time vector to the numeric reference date. 
% Important: As we're using numeric dates, the temporal resolution 
% does NOT play a role here! 
if strcmp(tme_unit{1}, 'hours')
    rel_dtes = double(time_in)/24;
elseif strcmp(tme_unit{1}, 'days')
    rel_dtes = double(time_in);
end

num_out = ref_dte_num + rel_dtes;
tme_out = datevec(num_out);

