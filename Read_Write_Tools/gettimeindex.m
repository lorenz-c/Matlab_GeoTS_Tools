function [start_indx, end_indx, count] = gettimeindex(fnme, start_date, end_date)
% The function opens the netcdf-file fnme and loads the time-dimension. The
% relative time axis is then transformed to absolute dates and the function
% looks for the indices of the start_date and end_date, which have both to
% be provided as [1 x 6] vectors, i.e. [yyyy mm dd hh mm ss]. The output of
% the function are the indices of the start- and end-date as well as the
% number of time-steps between these two dates. The start_indx and count
% can then be used as additional inputs in the netcdf2datastruct-function
% to load only a specific time period.
%--------------------------------------------------------------------------
% Input (required):
% - fnme        String with the file-name of the netcdf-file
% - start_date  First date of the desired time period as [1 x 6]-vector
% - end_date    Last date of the desired time period as [1 x 6]-vector
%
% Output
% - start_indx  Index of the first date in the netcdf-file
% - end_indx    Index of the last date in the netcdf-file
% - count       Number of time-steps between the start- and end-date
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         May 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------


% Open the netcdf-file 
ncid = netcdf.open(fnme);

% Read the time data
time_id  = netcdf.inqVarID(ncid, 'time');
time_vec = netcdf.getVar(ncid, time_id);

% Get the unit of the time-vector
att_id    = netcdf.inqAttID(ncid, time_id, 'units');
time_unit = netcdf.getAtt(ncid, time_id, 'units');

% Transform the relative dates to absolute dates
[tme_out, num_out] = reldate2absdate(time_vec, time_unit);

start_indx = find(num_out == datenum(start_date));
end_indx   = find(num_out == datenum(end_date));

if isempty(start_indx)
    warning('Start date out of bounds. Use first available date!')
    start_indx = 1;
end

if isempty(end_indx)
    warning('End date out of bounds. Use last available date!')
    end_indx = find(num_out == num_out(end));
end

count = end_indx - start_indx + 1;

netcdf.close(ncid)