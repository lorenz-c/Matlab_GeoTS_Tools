function [ind, count] = getlatlonindex(fnme, bbox)
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

% Read the latitude data
lat_id  = netcdf.inqVarID(ncid, 'lat');
lats    = netcdf.getVar(ncid, lat_id);

% Read the latitude data
lon_id  = netcdf.inqVarID(ncid, 'lon');
lons    = netcdf.getVar(ncid, lon_id);

% Close the netcdf-file
netcdf.close(ncid)

% Get the indices of the four corners
d_lon_0   = abs(lons - bbox(1));
d_lon_1   = abs(lons - bbox(2));
d_lat_0   = abs(lats - bbox(3));
d_lat_1   = abs(lats - bbox(4));

[tmp, ind_lon_0] = min(d_lon_0);
[tmp, ind_lon_1] = min(d_lon_1);
[tmp, ind_lat_0] = min(d_lat_0);
[tmp, ind_lat_1] = min(d_lat_1);

% Check the ordering of the latitudes
if ind_lat_1 < ind_lat_0
    tmp       = ind_lat_1;
    ind_lat_1 = ind_lat_0;
    ind_lat_0 = tmp;
end

% Check the ordering of the longitudes
if ind_lon_1 < ind_lon_0
    tmp       = ind_lon_1;
    ind_lon_1 = ind_lon_0;
    ind_lon_0 = tmp;
end


lat_count = ind_lat_1 - ind_lat_0 + 1;
lon_count = ind_lon_1 - ind_lon_0 + 1;

ind   = [ind_lon_0, ind_lon_1, ind_lat_0, ind_lat_1];
count = [lon_count, lat_count];


