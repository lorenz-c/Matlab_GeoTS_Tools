function otpt = multiregionaverage(inpt, region_map_struct, varargin)
% The function computes time-series of area-weighted means over selected 
% areas. These areas are defined by an id_map (a matrix where connected
% regions have the same id) and an area_id vector (or scalar), which 
% contains all the areas over which the input fields should be
% aggregated
%--------------------------------------------------------------------------
% Input:    - inpt      Cell array or structure which contains the input 
%                       fields.                
%           - id_map    Map which defines the different areas
%           - area_id   Vector (or scalar) which contains the ids of the 
%                       desired areas
%
%           - method    The function computes either a weighted mean
%                       (wmean), regular mean (mean), weighted sum (sum), 
%                       or regular sum (sum) over the regions defined in 
%                       the id_map
%                       default: wmean
%           - weights   default: true
%           - gridarea  Method for calculating the size of the gridcells
%                       Can be set to regular (default), haversine, cos, or
%                       vincenty
%           - dim_order Defines the order of the dimensions. If set to 1
%                       (default), the time-series will be arranged as 
%                       [time X regions] (which is required for e.g. #
%                       OPeNDAP-access). If set to 2, the time-series will 
%                       be arranged as [regions X time] (e.g. for the NCSS)

% Output:   - otpt      Matrix or structure array, which contains the 
%                       aggregated values over the regions which were 
%                       defined by the id_map and the area_ids.
%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   July 2011
%--------------------------------------------------------------------------
% Uses: area_wghts.m, cell2catchmat.m, copy_datastruct.m,
% create_datastruct.m, DataConvert.m
%--------------------------------------------------------------------------
% Updates: - 25.01.2013 For-loops removed, switched to matrix-based 
%                       computation 
%          - 14.10.2015 Added support for new data structures, brushed up
%                       help text and code
%--------------------------------------------------------------------------

% Checking input arguments and setting default values
pp = inputParser;
pp.addRequired('inpt', @(x) (isstruct(x) | ischar(x)));
pp.addRequired('map_params', @(x) (isstruct(x) | ischar(x) | ismatrix(x)));

pp.addParamValue('region_ids', []);                            
pp.addParamValue('method', 'nanmean')
pp.addParamValue('areaweights', true)
pp.addParamValue('gridarea', 'regular')
pp.addParamValue('calc_cntr', true)
pp.addParamValue('copy_map', false)
pp.addParamValue('dim_order', 1)
pp.addParamValue('nanpixl', false)

pp.parse(inpt, region_map_struct, varargin{:})

region_ids  = pp.Results.region_ids;
method      = pp.Results.method;
areaweights = pp.Results.areaweights;
gridarea    = pp.Results.gridarea;
calc_cntr   = pp.Results.calc_cntr;
copy_map    = pp.Results.copy_map;
dim_order   = pp.Results.dim_order;
nanpixl     = pp.Results.nanpixl;

clear pp

% 0. Load the dataset
if ischar(inpt)
    inpt = netcdf2datastruct(inpt);
end

% 2. Get the area-ids
if isempty(region_ids)
    if isstruct(region_map_struct)
        region_ids = region_map_struct.Data.regions;
    else
        region_ids = unique(region_map_struct(:));
        region_ids(isnan(region_ids)) = [];
    end
end

% 3. id_map 
if isstruct(region_map_struct)
    if inpt.Data.lat == flipud(region_map_struct.Data.lat)
        warning('Flipping region mask')
        region_map_struct.Data.region_map = ...
                                 flipud(region_map_struct.Data.region_map);
        region_map_struct.Data.lat = flipud(region_map_struct.Data.lat);
    end
    region_map = region_map_struct.Data.region_map;
else
    region_map = region_map_struct;
end

% Remove the "fixed" variables
vars = fieldnames(inpt.Variables);
for i = 1:length(vars)
    isfixed(i) = isfixedvar(vars{i});
end
vars(isfixed == 1) = [];
% Check for gridded variables 
isgrid = isgridvar(inpt, vars);
% Remove all non-gridded variables
vars(isgrid == 0) = [];

% Check the different variables for missing values
for i = 1:length(vars)
    if isfield(inpt.Variables.(vars{i}), 'FillValue')
        mval(i) = inpt.Variables.(vars{i}).FillValue;
    elseif isfield(inpt.Variables.(vars{i}), 'missing_value')
        mval(i) = inpt.Variables.(vars{i}).missing_value;
    else
        mval(i) = NaN;
    end
end
    
% Number of regions and time-steps
nr_regions = length(region_ids);

% Check the input for a time-variable. If no such variable is found, it is
% assumed that the input contains only a single map
if isfield(inpt.Variables, 'time')  
    nts = size(inpt.Data.time, 1);
else
    nts = 1;
end

% Get latitude and longitude
if isfield(inpt.Data, 'lat') & isfield(inpt.Data, 'lon')
    lat = inpt.Data.lat;
    lon = inpt.Data.lon;
else
    error('Input does not have lat- and lon-values')
end

% Create a binary mask to remove all "unwanted" elements
bin_mask = zeros(size(region_map));

% Set the pixels of all aggregation regions to 1
for i = 1:nr_regions
    bin_mask(region_map == region_ids(i)) = 1; 
end

% Add a comment about the averaging mask
if isstruct(region_map_struct)
    new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
         ': Spatial ', method, ' over ', region_map_struct.DataInfo.title];   
else
    new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                               ': Spatial ', method, ' over each region.'];    
end

% Go through the input dataset and set all missing elements to zero in the
% binary and the id map
for i = 1:length(vars)    
    % Create a single bin_mask for each time-step
    bin_mask_dta = repmat(bin_mask, 1, 1, nts);
    % Change the dimension order of the bin_mask (Time - lat - lon)
    bin_mask_dta = permute(bin_mask_dta, [3 1 2]);
    
    % Set the corresponding bin_mask-pixels to zero if the input-data
    % contains missing values
    if strcmp(method(1:3), 'nan')
        if isnan(mval(i))
            bin_mask_dta(isnan(inpt.Data.(vars{i})))  = 0;
        else
            bin_mask_dta(inpt.Data.(vars{i}) == mval) = 0;
        end
    end
    
    if nanpixl == true
        % Remove all pixels which contain at least one missing value during the
        % considered time-series
        bin_mask_dta                     = squeeze(sum(bin_mask_dta, 1));
        bin_mask_dta(bin_mask_dta < nts) = 0;
        bin_mask_dta(bin_mask_dta > 0)   = 1;
    else
        bin_mask_dta                     = squeeze(sum(bin_mask_dta, 1));
        bin_mask_dta(bin_mask_dta > 0)   = 1;
    end
    
    
    % Set all missing elements and unwanted areas to zero
    id_map_dta = region_map.*bin_mask_dta;

    % Create vectors from the masks
    bin_mask_dta_vec = bin_mask_dta(:);
    id_map_dta_vec   = id_map_dta(:);
    
    % Get the indices of valid elements
    cindx = find(bin_mask_dta_vec == 1);
    
    % Re-arrange the maps to vectors which contain only the non-zero 
    % elements of bin_mask
    id_map_dta_vec  = id_map_dta_vec(cindx);
    
    % For each region, a row in the matrix H is created. At this stage, H is
    % binary and defines if an element in the data matrix is located in the
    % current catchment or not.
    H = zeros(length(id_map_dta_vec), nr_regions);
    
    for j = 1:nr_regions
        tmp                                  = ones(length(id_map_dta_vec), 1);
        tmp(id_map_dta_vec ~= region_ids(j)) = 0;
        H(:, j)                              = tmp;   
    end
    
    % For weighted computations, a matrix A_mer is created which contains 
    % the areas of the pixels. This is used, depending on the chosen 
    % method, to apply area weights to the elements in the data matrix. As 
    % these are all linear operations (y=A*H), the weights are added to the 
    % H-matrix.
    if areaweights == true
        wghts = area_wghts(lat', lon, 'mat', gridarea);
        wghts = wghts(cindx);
    else
        wghts = ones(size(cindx));
    end
    
    for j = 1:nts
        tmp            = squeeze(inpt.Data.(vars{i})(j, :, :));
        tmp            = tmp(:);
        Data_mat(j, :) = tmp(cindx);
    end

    % Finally, compute the area average (sum, rms, variance, ...)
    otpt.Data.(vars{i}) = agg_data(Data_mat, H, wghts, method);  
    
    % Copy the variable meta-data from the input. 
    otpt.Variables.(vars{i})             = inpt.Variables.(vars{i});
    
    if dim_order == 1
        otpt.Variables.(vars{i}).dimensions  = {'time', 'regions'}; 
    elseif dim_order == 2
        otpt.Variables.(vars{i}).dimensions  = {'regions', 'time'};
        otpt.Data.(vars{i})                  = otpt.Data.(vars{i})';
    end
    
    % Add some short description of the calculations to the
    % source-attribute of each variable
    if isfield(otpt.Variables.(vars{i}), 'source')
        otpt.Variables.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                         otpt.Variables.(vars{i}).source]);
    else
        otpt.Variables.(vars{i}).source = new_srce;
    end
    
end

% Copy the meta-data from the gridded input-data
otpt.DataInfo = inpt.DataInfo;

% Copy the time-data from the inpt-data
if isfield(inpt.Variables, 'time')
    otpt.Dimensions.time = Inf;
    otpt.Variables.time  = inpt.Variables.time;
    otpt.Data.time       = inpt.Data.time;
    if isfield(inpt, 'TimeStamp')
        otpt.TimeStamp = inpt.TimeStamp;
    else
        otpt.TimeStamp = datenum(otpt.Data.time);
    end
end

% Copy the "region" variable from the region map to the output
if isstruct(region_map_struct)
    otpt = copyvars(otpt, region_map_struct, {'regions'});
    
    if isfield(region_map_struct.Variables, 'region_area')
        otpt = copyvars(otpt, region_map_struct, {'region_area'});
    end
    
    if isfield(region_map_struct.Variables, 'region_names')
        otpt = copyvars(otpt, region_map_struct, {'region_names'});
        otpt.Variables.region_names.nctype = 'NC_CHAR';
        % Check if the "max_str_length"-dimension is present                          
        if isfield(region_map_struct.Dimensions, 'max_str_length')
            otpt.Dimensions.max_str_length = ...
                               region_map_struct.Dimensions.max_str_length;
        else
            for i = 1:nr_regions
                mx_lngth(i) = ...
                            length(region_map_struct.Data.region_names{i});
            end
            otpt.Dimensions.max_str_length = max(mx_lngth);
        end
    end
    
    % Add the region_names to the coordinates-attribute of the variables
    for i = 1:length(vars)
        if isfield(otpt.Variables.(vars{i}), 'coordinates')
            otpt.Variables.(vars{i}).coordinates = ...
             strcat(otpt.Variables.(vars{i}).coordinates, ' region_names');
        else
            otpt.Variables.(vars{i}).coordinates = 'region_names';
        end
    end
    
else
    % If region_map_structure was only a matrix, we have to add the region
    % data and meta-data manually
    otpt.Dimensions.regions           = length(region_ids);
    otpt.Data.regions                 = region_ids;
    otpt.Variables.regions.long_name  = 'ID Number for each region';
    otpt.Variables.regions.short_name = 'station_id';
    otpt.Variables.regions.cf_role    = 'timeseries_id';
    otpt.Variables.regions.dimensions = {'regions'};
end


if copy_map == true
    % Create two dimensions for the index map
    otpt.Dimensions.lat_map = length(inpt.Data.lat);
    otpt.Dimensions.lon_map = length(inpt.Data.lon);
    
    % Copy the corresponding region map variable, latitude and longitude to 
    % the output
    if isstruct(region_map_struct)
        otpt = copyvars(otpt, region_map_struct, {'region_map'});
        % Latitude and longitude have to be renamed to lat_map and lon_map, 
        % as lat and lon might be used by the centers of the different areas
        otpt.Variables.lat_map = region_map_struct.Variables.lat;
        otpt.Variables.lon_map = region_map_struct.Variables.lon;
        otpt.Data.lat_map      = region_map_struct.Data.lat;
        otpt.Data.lon_map      = region_map_struct.Data.lon;
    else
        % If region_map_structure was only a matrix, we have to add some
        % meta-data manually.
        otpt.Data.region_map = region_map;
        otpt.Data.lat_map    = inpt.Data.lat;
        otpt.Data.lon_map    = inpt.Data.lon;
        
        otpt.Variables.region_map.long_name = ...
                                  'Region map with identification numbers';
        otpt.Variables.region_map.FillValue           = NaN;
        otpt.Variables.region_map.standard_name       = 'ID_Map';
        otpt.Variables.region_map.ancillary_variables = 'regions';
        otpt.Variables.lat_map                        = inpt.Variables.lat;
        otpt.Variables.lon_map                        = inpt.Variables.lon;
        otpt.Variables.lat_map.dimensions             = {'lat_map'};
        otpt.Variables.lon_map.dimensions             = {'lon_map'};                             
    end
        
    otpt.Variables.region_map.dimensions = {'lat_map', 'lon_map'};
end

% Calculate the area-centers
if calc_cntr == true
	if isvector(inpt.Data.lat) & isvector(inpt.Data.lon)
        % Number of lats and lons
        nr_lat = inpt.Dimensions.lat;
        nr_lon = inpt.Dimensions.lon;
        % Create two matrices with latitude and longitude values
        if isrow(inpt.Data.lat)
            lat_mat = repmat(inpt.Data.lat', nr_lon);
        else
            lat_mat = repmat(inpt.Data.lat, nr_lon);
        end
    
        if iscolumn(inpt.Data.lon)
            lon_mat = repmat(inpt.Data.lon', nr_lat);
        else
            lon_mat = repmat(inpt.Data.lon, nr_lat);
        end
    elseif ismatrix(inpt.Data.lat) & ismatrix(inpt.Data.lon)
        lat_mat = inpt.Data.lat;
        lon_mat = inpt.Data.lon;
    end
        
    % Loop over the regions
    for i = 1:nr_regions
        tmp_lat = lat_mat(region_map == region_ids(i));
        tmp_lon = lon_mat(region_map == region_ids(i));
        
        lat_cntr(i) = nanmean(tmp_lat);
        lon_cntr(i) = nanmean(tmp_lon);
    end
    
    % Add the required meta-data to the output
    otpt.Variables.lat.long_name  = 'Latitude of the center of the indexed areas';
    otpt.Variables.lat.units      = 'degrees_north';
    otpt.Variables.lat.dimensions = {'regions'};
    otpt.Data.lat                 = lat_cntr;

    otpt.Variables.lon.long_name  = 'Longitude of the center of the indexed areas';
    otpt.Variables.lon.units      = 'degrees_east';
    otpt.Variables.lon.dimensions = {'regions'};
    otpt.Data.lon                 = lon_cntr;    
    
    % Give each variable the lat- and lon-coordinates (Note: this might 
    % lead to some issues with external programs like, e.g., Panoply...) 
%     for i = 1:length(vars)
%         if isfield(otpt.Variables.(vars{i}), 'coordinates')
%             otpt.Variables.(vars{i}).coordinates = ...
%                   strcat(otpt.Variables.(vars{i}).coordinates, ' lat lon');
%         else  
%             otpt.Variables.(vars{i}).coordinates = 'lat lon';
%         end
%     end
end

% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                       '; MATLAB TS-Tools: spataverage.m'];
if isfield(inpt.DataInfo, 'history')           
    otpt.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                   inpt.DataInfo.history]);
else
    otpt.DataInfo.history = new_hist;
end


% Make sure that the output has the correct ordering
otpt = orderfields(otpt, {'DataInfo', ...
                          'Dimensions', ...
                          'Variables', ...
                          'Data', ...
                          'TimeStamp'});

                  
end    
        

















function Data_out = agg_data(dat_mtrx, mult, wghts, method)
    if strcmp(method, 'mean') | strcmp(method, 'nanmean')
        mult_wghted = bsxfun(@times, mult, wghts);
        area_sum    = sum(mult_wghted, 1);
        Data_out    = bsxfun(@rdivide, dat_mtrx*mult_wghted, area_sum);
    elseif strcmp(method, 'sum') | strcmp(method, 'nansum')
        mult_wghted = bsxfun(@times, mult, wghts);
        Data_out    = dat_mtrx*mult_wghted;
    elseif strcmp(method, 'sum_squared') | strcmp(method, 'nansum_squared')
        mult_wghted = bsxfun(@times, mult, wghts);
        Data_out    = dat_mtrx.^2*mult_wghted;
    elseif strcmp(method, 'rms') | strcmp(method, 'nanrms')
        mult_wghted = bsxfun(@times, mult, wghts);
        area_sum    = sum(mult_wghted, 1);
        tmp         = bsxfun(@rdivide, (dat_mtrx.^2)*mult_wghted, ...
                                                                 area_sum);
        Data_out    = sqrt(tmp);
    elseif strcmp(method, 'variance') | strcmp(method, 'nanvariance')
        % Compute the mean
        mult_wghted = bsxfun(@times, mult, wghts);
        V1          = sum(mult_wghted, 1);
        mn_dat      = bsxfun(@rdivide, dat_mtrx*mult_wghted, V1);
        % Create a matrix where each pixel within a region contains the
        % corresponding area average
        mn_pxl      = mn_dat*mult_wghted';
        % Remove the mean from the data and square the residuals
        inner       = (dat_mtrx - mn_pxl).^2;
        % Compute the sum by multiplicating the data matrix with the 
        % weights and the mask
        prod        = inner*mult_wghted;
        % Compute the squared sum of the weights
        V2 = sum(mult_wghted.^2);
        keyboard
        % Finally, compute the weighted variance
        Data_out = bsxfun(@times, prod, V1./(V1.^2 - V2)); 
    end 
end        
        
       
    
    
        
    


