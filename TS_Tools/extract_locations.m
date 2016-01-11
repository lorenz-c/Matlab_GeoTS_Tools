function out = extract_locations(inpt, sttn_dta, method, nghbrs)

if nargin < 4, nghbrs = 0; end
if nargin < 3, method = 'simple'; end

% Get the variables of inpt which are not fixed variables (dimensions, ...)
vars = fieldnames(inpt.Variables);

for i = 1:length(vars)
    isvar(i) = isfixedvar(vars{i});
end

% Remove all fixed variables from the variables-list
vars(isvar == 1) = [];

% Check the size of the lat- and lon-vectors
sze_lat = size(inpt.Data.lat);
if sze_lat(1) == 1 && sze_lat(2) > 1
    inpt.Data.lat = inpt.Data.lat';
end

sze_lon = size(inpt.Data.lon);
if sze_lon(1) == 1 && sze_lon(2) > 1
    inpt.Data.lon = inpt.Data.lon';
end

if isfield(inpt, 'TimeStamp')
    nts = length(inpt.TimeStamp);
else
    nts = size(inpt.Data.time, 1);
end

% First, compute the vertices of the input grids
if isfield(inpt.DataInfo, 'geospatial_lat_resolution') && ...
    isfield(inpt.DataInfo, 'geospatial_lon_resolution')

    % If the "meta-data" of intp contains the fields geospatial_lat_- and 
    % geospatial_lon_resolution, the function uses the attached values.
    lat_grid_top    = inpt.Data.lat + inpt.DataInfo.geospatial_lat_resolution/2;
    lat_grid_bottom = inpt.Data.lat - inpt.DataInfo.geospatial_lat_resolution/2;

    lon_grid_left   = inpt.Data.lon - inpt.DataInfo.geospatial_lon_resolution/2;
    lon_grid_right  = inpt.Data.lon + inpt.DataInfo.geospatial_lon_resolution/2;
    
else
    
    % If not, the resolution is computed from the data. Therefore, the inpt
    % variable must contain both lat and lon data.
    d_lat           = abs(inpt.Data.lat(1:end-1) - inpt.Data.lat(2:end));
    lat_grid_top    = inpt.Data.lat + 1/2*[d_lat(1); d_lat];
    lat_grid_bottom = inpt.Data.lat - 1/2*[d_lat;    d_lat(end)];
    
    d_lon           = abs(inpt.Data.lon(1:end-1) - inpt.Data.lon(2:end));
    lon_grid_left   = inpt.Data.lon - 1/2*[d_lon(1); d_lon];
    lon_grid_right  = inpt.Data.lon + 1/2*[d_lon;    d_lon(end)];
    
end


if strcmp(method, 'simple')
    
    for i = 1:length(sttn_dta.Data.lat)
   
        % Get the indices of the grid-cells which contain the given
        % locations. 
        lat_indx = find(sttn_dta.Data.lat(i) <= lat_grid_top  & ...
                        sttn_dta.Data.lat(i) >= lat_grid_bottom);
        lon_indx = find(sttn_dta.Data.lon(i) >= lon_grid_left & ...
                        sttn_dta.Data.lon(i) <= lon_grid_right);
                    
        if ~isempty(lat_indx) && ~isempty(lon_indx)
            lat_indx = [lat_indx(1) - nghbrs : 1 : lat_indx(end) + nghbrs];
            lon_indx = [lon_indx(1) - nghbrs : 1 : lon_indx(end) + nghbrs];
        
            % Check if any of the indices are beyond the map borders
            lat_above = find(lat_indx <= 0);
    
            if ~isempty(lat_above)
                lat_indx(lat_above) = lat_indx(lat_above(end) + 1);
            end
        
            lat_below = find(lat_indx > length(lat_grid_top));
        
            if ~isempty(lat_below)
                lat_indx(lat_below) = lat_indx(lat_below(1) - 1);
            end
        
            lon_left = find(lon_indx <= 0);
        
            if ~isempty(lon_left)
                if min(lon_left) == -180
                    lon_indx(lon_left) = lon_indx(lon_left) + length(lon_grid_right);
                else
                    lon_indx(lon_left) = lon_indx(lon_left(end) + 1);
                end
            end
        
            lon_right = find(lon_indx > length(lon_grid_right));
        
            if ~isempty(lon_right)
                if max(lon_right) == 180
                    lon_indx(lon_right) = lon_indx(lon_right) - 360/d_lon(1);
                else
                    lon_indx(lon_right) = lon_indx(lon_right(1) - 1);
                end
            end
        end
            
        for j = 1:length(vars)
            if ~isempty(lat_indx) && ~isempty(lon_indx)
                % Extract the data using the previously calculated indices
                tmp_data = inpt.Data.(vars{j})(:, lat_indx, lon_indx);
    
                % Compute the mean over the sub-matrices
                data_out = nanmean(nanmean(tmp_data, 3), 2);
    
                out.Data.(vars{j})(:, i) = data_out;    
            else
                out.Data.(vars{j})(1:nts, i) = NaN;
            end
        end
    end
    
    
elseif strcmp(method, 'IDW')

    for i = 1:length(sttn_dta.Data.lat)
    
        % Get the indices of the grid-cells which contain the given
        % locations. 
        lat_indx = find(sttn_dta.Data.lat(i) <= lat_grid_top  & sttn_dta.Data.lat(i) >= lat_grid_bottom);
        lon_indx = find(sttn_dta.Data.lon(i) >= lon_grid_left & sttn_dta.Data.lon(i) <= lon_grid_right);
        
        if ~isempty(lat_indx) & ~isempty(lon_indx)
            lat_indx = [lat_indx(1) - nghbrs : 1 : lat_indx(end) + nghbrs];
            lon_indx = [lon_indx(1) - nghbrs : 1 : lon_indx(end) + nghbrs];
            
            % Check if any of the indices are beyond the map borders
            lat_above = find(lat_indx <= 0);
    
            if ~isempty(lat_above)
                lat_indx(lat_above) = lat_indx(lat_above(end) + 1);
            end
        
            lat_below = find(lat_indx > length(lat_grid_top));
        
            if ~isempty(lat_below)
                lat_indx(lat_below) = lat_indx(lat_below(1) - 1);
            end
        
            lon_left = find(lon_indx <= 0);
        
            if ~isempty(lon_left)
                lon_indx(lon_left) = lon_indx(lon_left) + length(lon_grid_right);
            end
        
            lon_right = find(lon_indx > length(lon_grid_right));
        
            if ~isempty(lon_right)
                lon_indx(lon_right) = lon_indx(lon_right) - length(lon_grid_right);
            end
             
            % Create two matrices which contain the lats and lons of the input
            % grid
            tmp_lat  = repmat(inpt.Data.lat(lat_indx), 1, length(lon_indx));
            tmp_lon  = repmat(inpt.Data.lon(lon_indx)', length(lat_indx), 1);
    
            % Compute the distance between the grid cells and the locations of
            % the station data
            dist = haversine(tmp_lat, tmp_lon, ...
                        sttn_dta.Data.lat(i), sttn_dta.Data.lon(i), ...
                        6371.137, false, false);
    
            % Add "third" dimension to the dist-matrix
            dist = reshape(dist, [1 size(dist)]);
            dist = repmat(dist, nts, 1, 1);
        
            for j = 1:length(vars)
                % Extract the data using the previously calculated indices
                tmp_data = inpt.Data.(vars{j})(:, lat_indx, lon_indx);
            
                % Remove the missing elements in the data from the dist-array 
                dist(isnan(tmp_data)) = NaN;
            
                % Multiply the data with the distances
                data_wghted = tmp_data.*dist;
            
                % Compute the final weighted mean over all grid-cells
                data_out = nansum(nansum(data_wghted, 3), 2)./ ...
                       nansum(nansum(dist, 3), 2);
    
                out.Data.(vars{j})(:, i) = data_out;    
            end
        else
            for j = 1:length(vars)
                out.Data.(vars{j})(:, i) = NaN(nts, 1);
            end
        end           
    end
end

% Copy the meta-data from the input
out.DataInfo  = inpt.DataInfo;

for i = 1:length(vars)
    % Copy the variable attributes from the input
    out.Variables.(vars{i})             = inpt.Variables.(vars{i});
    % ...but add the new dimensions
    out.Variables.(vars{i}).dimensions  = {'time', 'stations'};
    % ...and coordinates
    if isfield(sttn_dta.Variables, 'alt')
        out.Variables.(vars{i}).coordinates = 'station_name lat lon alt';
    else
        out.Variables.(vars{i}).coordinates = 'station_name lat lon';
    end
end

% Use the time-data from the input variable
out = copyvars(out, inpt, {'time'});
out.TimeStamp = datenum(out.Data.time);

% Replace the 3D-Grid dimensions with the station dimensions from the
% sttn_dta variable
out.Dimensions = sttn_dta.Dimensions;

% Copy the station-locations, etc. to the output variable
out = copyvars(out, sttn_dta, {'lat', 'lon', 'stations', 'station_name', 'alt'});

% Update the file history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
    '; MATLAB TS-Tools: extract_locations.m; ', method, ...
    ' average using ', num2str(nghbrs), ' neighbours.'];    

if isfield(out.DataInfo, 'history')
    out.DataInfo.history = sprintf([new_hist, ' \n', out.DataInfo.history]);
else
    out.DataInfo.history = new_hist;
end

             




