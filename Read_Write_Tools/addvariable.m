function otpt = addvariable(inpt, varnames, vartype, varprec)
% The function creates an empty datastructure with the respective fields
% for the data, variables, and metadata. Depending on the chosen variable
% type, different attributes, dimensions, and empty data arrays are added
% to the variable.
%--------------------------------------------------------------------------
% INPUT:
% - varnames    Cell-array where each cell contains a variable name. For
%               each variable name, a separate variable and empty data 
%               array created with the respective dimensions, attributes,
%               etc.
% - vartype     Type of the variable. Possible choices are '1d_stations',
%               '2d_stations', '1d_regions', '2d_regions', '2d_grids', 
%               '3d_grids', and '4d_grids'
% - varprec     String which defines the prescission of the variable.
%               Possible choices are 'double', 'single', 'int32', 'int16',
%               and 'int8'. More will follow.
%--------------------------------------------------------------------------
% OUTPUT:
% - out         Matlab datastructure
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: addvariable.m
%--------------------------------------------------------------------------


otpt = inpt;
if nargin < 4, varprec = 'double'; end

if isstruct(inpt.Variables)
    vars = fieldnames(inpt.Variables);
else
    vars = {''};
end

% 1. Check if lat, lon, time are already present in inpt
if ismember(vartype, {'2d_grids', '3d_grids', '4d_grids'})
    
    % Check for latitude
    if ~ismember(vars, 'lat')

        % Add latitude to inpt
        otpt.Variables.lat.long_name     = 'Latitude';
        otpt.Variables.lat.standard_name = 'latitude';
        otpt.Variables.lat.units         = 'degrees_north';
        otpt.Variables.lat.dimensions    = {'lat'};
        
        otpt.Dimensions.lat = double.empty(0);
        otpt.Data.lat       = double.empty(0);
    end
    
    % Check for longitude
    if ~ismember(vars, 'lon')
        % Add longitude to inpt
        otpt.Variables.lon.long_name     = 'Longitude';
        otpt.Variables.lon.standard_name = 'latitude';
        otpt.Variables.lon.units         = 'degrees_east';
        otpt.Variables.lon.dimensions    = {'lon'};
        
        otpt.Dimensions.lon = double.empty(0);
        otpt.Data.lon       = double.empty(0);
    end
    
elseif strcmp(vartype, '2d_stations')
    
    % Check for stations
    if ~ismember(vars, 'stations')
        otpt.Variables.stations.long_name     = char.empty(0);
        otpt.Variables.stations.standard_name = 'station_id';
        otpt.Variables.stations.cf_role       = 'timeseries_id';
        otpt.Variables.stations.ancillary_variables = char.empty(0);
        otpt.Variables.stations.dimensions    = {'stations'};
        
        otpt.Dimensions.stations = single.empty(0);  
        otpt.Data.stations       = int64.empty(0);
    end
        
elseif strcmp(vartype, '2d_regions')
    
    % Check for stations
    if ~strcmp(vars, 'regions')
        otpt.Variables.regions.long_name     = char.empty(0);
        otpt.Variables.regions.standard_name = 'station_id';
        otpt.Variables.regions.cf_role       = 'timeseries_id';
        otpt.Variables.regions.ancillary_variables = char.empty(0);
        otpt.Variables.regions.dimensions    = {'stations'};
        
        otpt.Dimensions.regions = single.empty(0);  
        otpt.Data.regions       = int64.empty(0);
    end
      
end



if ismember(vartype, {'1d_stations', '1d_regions'})
    gridsize = [0];
elseif ismember(vartype, {'2d_grids', '2d_stations', '2d_regions'}) 
    gridsize = [0 0];
elseif ismember(vartype, '3d_grids')
    gridsize = [0 0 0];
elseif ismember(vartype, '4d_grids')
    gridsize = [0 0 0 0];
end
 

if ismember(vartype, {'1d_stations', '2d_stations', '1d_regions', ...
                                     '2d_regions', '3d_grids', '4d_grids'})
        
    if ~ismember(vars, 'time')
        % Add latitude to inpt
        otpt.Variables.time.long_name     = 'Time';
        otpt.Variables.time.standard_name = 'time';
        otpt.Variables.time.units         = 'yyyy-MM-dd HH:mm:ss';
        otpt.Variables.time.calendar      = 'standard';
        otpt.Variables.time.dimensions    = {'time'};
        
        otpt.Dimensions.time = double.empty(0);
        otpt.Data.time       = double.empty(0);
        otpt.TimeStamp       = double.empty(0);
    end    
end
    
if strcmp(vartype, '4d_grids')
    if ~strcmp(vars, 'levels')
        % Add latitude to inpt
        otpt.Variables.time.long_name     = char.empty(0);
        otpt.Variables.time.standard_name = 'levels';
        otpt.Variables.time.units         = char.empty(0);
        otpt.Variables.time.dimensions    = {'levels'};
        
        otpt.Dimensions.levels = [];
        otpt.Data.levels       = [];

    end    
end

if ismember(vartype, {'1d_stations', '1d_regions'})
    datadims = {'time'};  
    otpt.DataInfo.cdm_data_type = 'Station';
    otpt.DataInfo.featureType   = 'timeSeries';  
elseif strcmp(vartype, '2d_stations')
    datadims = {'stations', 'time'};
    otpt.DataInfo.cdm_data_type = 'Station';
    otpt.DataInfo.featureType   = 'timeSeries';  
elseif strcmp(vartype, '2d_regions')
    datadims = {'regions', 'time'};
    otpt.DataInfo.cdm_data_type = 'Station';
    otpt.DataInfo.featureType   = 'timeSeries';
elseif strcmp(vartype, '2d_grids')
    datadims = {'lat', 'lon'};
    otpt.DataInfo.cdm_data_type = 'Grid';
elseif strcmp(vartype, '3d_grids')
    datadims = {'time', 'lat', 'lon'};
    otpt.DataInfo.cdm_data_type = 'Grid';
elseif strcmp(vartype, '4d_grids')
    datadims = {'time', 'levels', 'lat', 'lon'};
    otpt.DataInfo.cdm_data_type = 'Grid';
end

for i = 1:length(varnames)   
    tmp = struct('long_name', char.empty(0), ...
                 'standard_name', char.empty(0), ...
                 'units', char.empty(0), ...
                 'FillValue', NaN, ...
                 'dimensions', {datadims});
             
    otpt.Variables.(varnames{i}) = tmp;
    
    if strcmp(varprec, 'double')
        otpt.Data.(varnames{i}) = double.empty(gridsize);
    elseif strcmp(varprec, 'single')
        otpt.Data.(varnames{i}) = single.empty(gridsize);
        otpt.Variables.(varnames{i}).nctype = 'NC_FLOAT';
    elseif strcmp(varprec, 'int32')
        otpt.Data.(varnames{i}) = int32.empty(gridsize);
        otpt.Variables.(varnames{i}).nctype = 'NC_INT';
    elseif strcmp(varprec, 'int16')
        otpt.Data.(varnames{i}) = int16.empty(gridsize);
        otpt.Variables.(varnames{i}).nctype = 'NC_SHORT';
    elseif strcmp(varprec, 'int8')
        otpt.Data.(varnames{i}) = int8.empty(gridsize);
        otpt.Variables.(varnames{i}).nctype = 'NC_BYTE';    
    end
end
    
    
    
    
    
    
    
    
        
    