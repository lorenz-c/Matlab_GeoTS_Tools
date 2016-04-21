function station_struct = indomain(station_struct, grid_lat, grid_lon, mask_dta)
% The function selects the stations in the input variable station_struct
% which are located within a specific domain. At first, the function
% computes a bounding box, depending on the coordinates in grid_lat and
% grid_lon. Only the stations within this box are added to the output data.
% If furthermore a mask (mask_dta) is provided (which contains real values
% in the region of interest and NaNs elsewhere), the function selects only
% the stations which are located within this region.
%--------------------------------------------------------------------------
% INPUT:
% - station_struct    Datastructure with station data
% - grid_lat,         Latitude and longitude of the region of interest.
%   grid_lon          Note that grid_lat and grid_lon are not the 
%                     coordinates of a polygon, but the coordinates of grid 
%                     cell centers.
% - mask_dta          Matrix with the dimension 
%                     [length(grid_lat) x length(grid_lon)]. The region of
%                     interest is identified by real values. The grid cells
%                     outside of this region must be set to NaN.
%--------------------------------------------------------------------------
% OUTPUT:
% - station_struct    Datastructure with the stations, which are located in
%                     the region of interest.
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: isstationvar.m, istimevar.m, getdimpos.m
%--------------------------------------------------------------------------     
if nargin < 4, mask_dta = []; end

% Get lats and lons of the station data
stations_lat = station_struct.Data.lat;
stations_lon = station_struct.Data.lon;

% Get all "station" variables (i.e. vars with the "station"-dimension)
vars                 = fieldnames(station_struct.Variables);
isstation            = isstationvar(station_struct, vars);
vars(isstation == 0) = [];

% Get all stations which are located within the bounding box
station_indx = find(stations_lat <= max(grid_lat) & ...
                    stations_lat >= min(grid_lat) & ...
                    stations_lon >= min(grid_lon) & ...
                    stations_lon <= max(grid_lon));
               
% Update the size of the stations-dimension
station_struct.Dimensions.stations = length(station_indx);

% Get all variables with a time dimension                
tme_vars     = istimevar(station_struct, vars); 

% Get the "position" of the time dimension
dimpos = getdimpos(station_struct, vars, 'time');

% Loop over the station variables
for i = 1:length(vars)
    % Check for the time dimension
    if tme_vars(i) == 1
        % Check if "time" is the first or second dimension of the current
        % variable
        if dimpos(i) == 1
            station_struct.Data.(vars{i}) = ...
                            station_struct.Data.(vars{i})(:, station_indx);
        elseif dimpos(i) == 2
            station_struct.Data.(vars{i}) = ...
                            station_struct.Data.(vars{i})(station_indx, :);
        end
    else
        % Get the data with the indices "station_indx"
        station_struct.Data.(vars{i}) = ...
                               station_struct.Data.(vars{i})(station_indx);
    end
end

if ~isempty(mask_dta)
    % Generate a binary mask 
    mask                  = ones(size(mask_dta));
    mask(isnan(mask_dta)) = 0;
    
    % Compute the coordinates of the left, right, bottom, and top edges of
    % the grid cells
    d_lat           = abs(grid_lat(1:end-1) - grid_lat(2:end));
    lat_grid_top    = grid_lat + 1/2*[d_lat(1); d_lat];
    lat_grid_bottom = grid_lat - 1/2*[d_lat;    d_lat(end)];
    
    d_lon           = abs(grid_lon(1:end-1) - grid_lon(2:end));
    lon_grid_left   = grid_lon - 1/2*[d_lon(1); d_lon];
    lon_grid_right  = grid_lon + 1/2*[d_lon;    d_lon(end)];
    
    % Generate a grid from the pixel edges with the same size as the binary 
    % mask
    [left_lon,     top_lat] = meshgrid(lon_grid_left, lat_grid_top);
    [right_lon, bottom_lat] = meshgrid(lon_grid_right, lat_grid_bottom);
    
    % Get the new station coordinates
    stations_lat = station_struct.Data.lat;
    stations_lon = station_struct.Data.lon;
    
    % Get the pixels which are located in the domain
    left_lon   = left_lon(mask == 1);
    right_lon  = right_lon(mask == 1);
    top_lat    = top_lat(mask == 1);
    bottom_lat = bottom_lat(mask == 1);
    
    % Create an empty station_indx vector
    station_indx = [];
    
    for i = 1:length(left_lon)
        % Check if there are one or more pixels within each pixel of the
        % domain.
        act_cell = find(stations_lat <= top_lat(i)    & ...
                        stations_lat >= bottom_lat(i) & ...
                        stations_lon >= left_lon(i)   & ...
                        stations_lon <= right_lon(i));
                    
        if ~isempty(act_cell)
            % Save the indices of the stations within the current pixels
            station_indx = [station_indx; act_cell];
        end
    end
    
    % As some stations might be located directly on an edge between two
    % pixels, the same index might appear twice in the station_indx.
    % Therefore, we have to remove multiple occurences of the same index. 
    station_indx = unique(station_indx);
    
    % Same loop as before, but with a new station_indx variable
    for i = 1:length(vars)
        if tme_vars(i) == 1
            if dimpos(i) == 1
                if strcmp(vars{i}, 'time')
                end
                station_struct.Data.(vars{i}) = ...
                            station_struct.Data.(vars{i})(:, station_indx);
            elseif dimpos(i) == 2
                station_struct.Data.(vars{i}) = ...
                            station_struct.Data.(vars{i})(station_indx, :);
            end
        else
            station_struct.Data.(vars{i}) = ...
            station_struct.Data.(vars{i})(station_indx);
        end
    end
end

% Compute the length of the station names and use that value for the
% max_str_length dimension
if isfield(station_struct.Variables, 'station_name')
    for i = 1:length(station_struct.Data.station_name)
        lnth(i) = length(station_struct.Data.station_name{i});
    end
    
    station_struct.Dimensions.max_str_length = max(lnth);
end
             
% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                    '; MATLAB TS-Tools: indomain.m'];
if isfield(station_struct.DataInfo, 'history')           
    station_struct.DataInfo.history = sprintf([new_hist, ' \n', ...
                                         station_struct.DataInfo.history]);
else
    station_struct.DataInfo.history = new_hist;
end             
             
             
        
    

