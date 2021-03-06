function otpt = spatialaverage(inpt, varargin)
% TBA
%--------------------------------------------------------------------------
% INPUT:
% - inpt        CF datastructure
% - vars        Cell-array of a list with variables, for which the averages
%               should be computed. If left empty, the function calculates
%               the averages for each gridded variable in inpt.
% - method      Defines the quantity of the output. Possible selections are
%               mean (default), sum, sum_squared, rms, variance
% - weights     Defines if spatial weights depending on the pixel areas 
%               should be applied.
% - gridarea    Defines the method for computing the pixel areas. Possible
%               selections are regular (default), haversine, cos, vincenty
%               See area_wghts.m for the description of the different
%               methods.
%--------------------------------------------------------------------------
% OUTPUT:
% - otpt        Datastructure with the spatial averages (or sums, ...) of
%               inpt.
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: area_wghts.m, isfixedvar.m, isgridvar.m
%--------------------------------------------------------------------------
% References
% - https://rulesofreason.wordpress.com/2012/02/13/weighted-variance-and-weighted-coefficient-of-variation/
%--------------------------------------------------------------------------


% Checking input arguments and setting default values
pp = inputParser;
pp.addRequired('inpt', @(x) (isstruct(x) | ischar(x)));

pp.addParamValue('vars', 'all');                  
pp.addParamValue('method', 'mean')
pp.addParamValue('areaweights', true)
pp.addParamValue('gridarea', 'regular')


pp.parse(inpt, varargin{:})

vars        = pp.Results.vars;
method      = pp.Results.method;
areaweights = pp.Results.areaweights;
gridarea    = pp.Results.gridarea;


clear pp

% 0. Load the dataset
if ischar(inpt)
    inpt = netcdf2datastruct(inpt);
end

% Get latitude and longitude
if isfield(inpt.Data, 'lat') & isfield(inpt.Data, 'lon')
    lat = inpt.Data.lat;
    lon = inpt.Data.lon;
else
    error('Input does not have lat- and lon-values')
end

if strcmp(vars, 'all')
    % Remove the "fixed" variables
    vars = fieldnames(inpt.Variables);
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
    end
    vars(isfixed == 1) = [];
    % Check for gridded variables
    isgrid = isgridvar(inpt, vars);
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



for i = 1:length(vars)
    % Check for missing values
    if isfield(inpt.Variables.(vars{i}), 'FillValue')
        mval = inpt.Variables.(vars{i}).FillValue;
    elseif isfield(inpt.Variables.(vars{i}), 'missing_value')
        mval = inpt.Variables.(vars{i}).missing_value;
    end
    
    % Create an empty region map
    region_map = zeros(size(inpt.Data.(vars{i})));
    
    % Set each grid cell with valid data to 1
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
    
    % Remove these pixels from the region map
    region_map = bsxfun(@times, region_map, nr_dta);

    if areaweights == true
        % Compute the area of the grid cells
        wghts = area_wghts(lat', lon, 'mat', gridarea);
        % Add a dummy dimension
        wghts = reshape(wghts, [1 size(wghts)]);    
        % Remove all unwanted elements from wghts
        wghts = wghts.*region_map(1, :, :);
    else
        % If no weights should be applied, we use simply the first region
        % map
        wghts = region_map(1, :, :);
    end
    
    
    % Compute the sum of the weights. If no weights should be applied,
    % area_sum is simply the number of grid cells with real data.
    area_sum = nansum(nansum(wghts, 3), 2);
    
    if strcmp(method, 'mean') 
        % Multiply the data with the weights
        fld      = bsxfun(@times, inpt.Data.(vars{i}), wghts);
        Data_out = nansum(nansum(fld, 3), 2)./area_sum;
    elseif strcmp(method, 'sum') 
        % Multiply the data with the weights
        fld      = bsxfun(@times, inpt.Data.(vars{i}), wghts);
        Data_out = nansum(nansum(fld, 3), 2);
    elseif strcmp(method, 'sum_squared') 
        % Multiply the data with the weights
        fld      = bsxfun(@times, inpt.Data.(vars{i}).^2, wghts);
        Data_out = nansum(nansum(fld, 3), 2);
    elseif strcmp(method, 'rms') 
        % Multiply the data with the weights
        fld      = bsxfun(@times, inpt.Data.(vars{i}).^2, wghts);
        Data_out = sqrt(nansum(nansum(fld, 3), 2)./area_sum);
    elseif strcmp(method, 'variance') 
        % Compute the sum of the weights
        V1 = nansum(nansum(wghts, 3), 2);
        V2 = nansum(nansum(wghts.^2, 3), 2);
        % Compute the mean
        mn = nansum(nansum(fld, 3), 2)./area_sum;
        % Compute the anomalies
        anom = bsxfun(@minus, fld, mn);
        % Multiply the squared anomalies with the weights
        anom = bsxfun(@times, anom.^2, wghts);
        % Compute the variance
        Data_out = (V1./(V1.^2 - V2)).*nansum(nansum(anom, 3), 2);
    end
    
    % Copy the averages to the output
    otpt.Data.(vars{i}) = Data_out;
    
    % Copy the variable meta-data from the input
    otpt.Variables.(vars{i})             = inpt.Variables.(vars{i});
    % Set the dimensions-attribute to "time", as this is the only dimension
    % of the new data.
    otpt.Variables.(vars{i}).dimensions  = {'time'}; 
    
    
    % Add the applied method to the source-attribute
    new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                        ': Spatial ', method];   
    
    if isfield(otpt.Data.(vars{i}), 'source')
        otpt.Variables.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                         otpt.Variables.(vars{i}).source]);
    else
        otpt.Variables.(vars{i}).source = new_srce;
    end
    
end    


% Copy the global meta-data from the gridded input-data
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
    otpt.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                   inpt.DataInfo.history]);
else
    otpt.DataInfo.history = new_hist;
end



        

        
       
    
    
        
    


