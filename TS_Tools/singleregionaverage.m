function otpt = spatialaverage(inpt, varargin)
% The function computes time-series of area-weighted means over selected 
% areas. These areas are defined by an id_map (a matrix where connected
% regions have the same id) and an area_id vector (or scalar), which 
% contains all the areas over which the input fields should be
% aggregated
%--------------------------------------------------------------------------
% Input:    - inpt      Cell array or structure which contains the input 
%                       fields.                
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
% References
% - https://rulesofreason.wordpress.com/2012/02/13/weighted-variance-and-weighted-coefficient-of-variation/

% Checking input arguments and setting default values
pp = inputParser;
pp.addRequired('inpt', @(x) (isstruct(x) | ischar(x)));

pp.addParamValue('vars', 'all');                  
pp.addParamValue('method', 'nanmean')
pp.addParamValue('areaweights', true)
pp.addParamValue('gridarea', 'regular')
pp.addParamValue('calc_cntr', true)
pp.addParamValue('copy_map', true)

pp.parse(inpt, varargin{:})

vars        = pp.Results.vars;
method      = pp.Results.method;
areaweights = pp.Results.areaweights;
gridarea    = pp.Results.gridarea;
calc_cntr   = pp.Results.calc_cntr;
copy_map    = pp.Results.copy_map;

clear pp

% 0. Load the dataset
if ischar(inpt)
    inpt = netcdf2datastruct(inpt);
end

if strcmp(vars, 'all')
    % Remove the "fixed" variables
    vars = fieldnames(inpt.Variables);

    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
    end
    vars(isfixed == 1) = [];

    % Check the remaining variables if they have a lat- and lon-dimension
    for i = 1:length(vars)
        isgrid(i) = find(ismember(inpt.Variables.(vars{i}).dimensions, ...
           'lat')) & find(ismember(inpt.Variables.(vars{i}).dimensions, ...
           'lon'));
    end

    % Remove variables which do not have the lat- and lon-dimension
    vars(isgrid == 0) = [];
end
  

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

for i = 1:length(vars)
    % Check for missing values
    if isfield(inpt.Variables.(vars{i}), 'FillValue')
        mval = inpt.Variables.(vars{i}).FillValue;
    elseif isfield(inpt.Variables.(vars{i}), 'missing_value')
        mval = inpt.Variables.(vars{i}).missing_value;
    end

    region_map = zeros(size(inpt.Data.(vars{i})));
    if isnan(mval)
        region_map(~isnan(inpt.Data.(vars{i}))) = 1;
    else
        region_map(inpt.Data.(vars{i}) ~= mval) = 1;
    end
    region_id = 1;
    
    % Compute the number of valid data for each pixel
    nr_dta = sum(region_map, 1);
    
    % Remove the pixels which do not have a "continuous" time series
    nr_dta(nr_dta < nts)  = NaN;
    nr_dta(nr_dta == nts) = 1;
    
    region_map = bsxfun(@times, region_map, nr_dta);
    
    % For weighted computations, a matrix A_mer is created which contains 
    % the areas of the pixels. This is used, depending on the chosen 
    % method, to apply area weights to the elements in the data matrix. As 
    % these are all linear operations (y=A*H), the weights are added to the 
    % H-matrix.
    if areaweights == true
        wghts = area_wghts(lat', lon, 'mat', gridarea);
        wghts = reshape(wghts, [1 length(lat) length(lon)]);        
        wghts = wghts.*region_map(1, :, :);
    else
        wghts = eye(1, length(lat), length(lon));
        wghts = wghts.*region_map(1, :, :);
    end
    
    fld      = bsxfun(@times, inpt.Data.(vars{i}), wghts);
    area_sum = nansum(nansum(wghts, 3), 2);
    
    if strcmp(method, 'mean') | strcmp(method, 'nanmean')
        Data_out = nansum(nansum(fld, 3), 2)./area_sum;
    elseif strcmp(method, 'sum') | strcmp(method, 'nansum')
        Data_out = nansum(nansum(fld, 3), 2);
    elseif strcmp(method, 'sum_squared') | strcmp(method, 'nansum_squared')
        Data_out = nansum(nansum(fld.^2, 3), 2);
    elseif strcmp(method, 'rms') | strcmp(method, 'nanrms')
        Data_out = nansum(nansum(fld.^2, 3), 2)./area_sum;
    elseif strcmp(method, 'variance') | strcmp(method, 'nanvariance')

        V1 = nansum(nansum(wghts, 3, 2));
        V2 = nansum(nansum(wghts.^2, 3, 2));
        
        % Compute the mean
        mn = nansum(nansum(fld, 3), 2)./area_sum;
        % Compute the anomalies
        anom = bsxfun(@minus, fld, mn);
        % Multiply the squared anomalies with the weights
        anom = bsxfun(@times, anom.^2, wghts);
        
        Data_out = (V1./(V1.^2 - V2)).*sum(sum(anom, 3, 2));
    end
    % Copy the averages to the output
    otpt.Data.(vars{i}) = Data_out;
    
    % Copy the variable meta-data from the input
    otpt.Variables.(vars{i})             = inpt.Variables.(vars{i});
    otpt.Variables.(vars{i}).dimensions  = {'time'}; 
    
    
    % Add a comment about the averaging mask
    new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                                     ': Spatial average '];   

    if isfield(otpt.Data.(vars{i}), 'source')
        otpt.Data.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                                    inpt.DataInfo.source]);
    else
        otpt.Data.(vars{i}).source = new_srce;
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

% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                    '; MATLAB TS-Tools: spatialaverage.m'];
if isfield(inpt.DataInfo, 'history')           
    otpt.DataInfo.history = sprintf([new_hist, ' \n', inpt.DataInfo.history]);
else
    otpt.DataInfo.history = new_hist;
end



        

        
       
    
    
        
    


