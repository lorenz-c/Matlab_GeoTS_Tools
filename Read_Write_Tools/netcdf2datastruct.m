function [otpt] = netcdf2datastruct(fnme, tme_trafo, nan_mval, vars, period, latlons)
% The function converts a netcdf-file into a datastructure. The function
% generally reads all variables, dimensions, attributes, etc. from the
% input file "fnme".
%--------------------------------------------------------------------------
% Input (required):
% - fnme        File name of the netCDF-file 
% - tme_trafo   Transform the time unit from, e.g., "days since ..." into
%               years, months, days, ...
%               Note: This works only if the time-variable in the
%               netCDF-file is given in units of "days since yyyy-mm-dd
%               hh:mm:ss" or "hours since yyyy-mm-dd hh:mm:ss"
%       tru        TBA: Support for time zones (e.g, UTC +01)
% - nan_mval    Boolean variable: if set to true, the function converts the
%               missing values into NaNs.
%               provided, the _FillValue in the netCDF-file is used as 
%               identifier. However, it is more convenient to work with 
%               matlab's NaNs.
% - vars        Cell-array of variables which should be loaded. The
%               function automatically checks for all mandatory dimensions 
%               and the corresponding variables. 
% - period      [1 x 2]-Vector which can be provided for loading only a 
%               certain period of time. The first element must be the index
%               of the first desired date and the second element must be
%               the number of time-steps, which should be loaded.
% - latlons     [1 x 4]-Vector which can be provided for loading only a
%               certain region. The first (second) element must be the
%               longitude (latitude) index of the lower left corner. The
%               third and fourth elements are the number of longitudes and
%               latitudes.
% Output:
% - otpt        Matlab datastructure
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
% Call: 
% [otpt] = netcdf2datastruct(fnme, tme_trafo, nan_mval, vars, period, latlons)
if nargin < 6, latlons = []; end

if nargin < 5, period = []; end

if nargin < 4, vars = 'all'; end

% Set all missing values to NaN
if nargin < 3, nan_mval = true; end

% The function transforms the time-vector by default
if nargin < 2, tme_trafo = true; end

% Open the netCDF-file
ncid = netcdf.open(fnme, 'NC_NOWRITE');

% Get general information about the netCDF-file
[numdims, numvars, numglobalatts, unlimdimid] = netcdf.inq(ncid);

%--------------------------------------------------------------------------
%                           GLOBAL ATTRIBUTES
%--------------------------------------------------------------------------
% Get the global ID
globID = netcdf.getConstant('NC_GLOBAL');

% Loop over the global attributes
for i = 1:numglobalatts
    % Get global attribute names of the corresponding IDs
    attname  = netcdf.inqAttName(ncid, globID, i-1);
    
    % Check for dashs and dots in the attribute name. As Matlab does not
    % allow these characters in the name of a structure, they have to be
    % replaced with something else.
    dash_pos  = find(ismember(attname, '-'));
    dot_pos   = find(ismember(attname, '.'));
    white_pos = find(ismember(attname, ' '));
    point_pos = find(ismember(attname, ':'));

    
    attname_mod = attname;
    
    if ~isempty(dash_pos)
        for j = 1:length(dash_pos)
            attname_mod(dash_pos(j)) = '_';
        end
    end
    
    if ~isempty(dot_pos)
        for j = 1:length(dot_pos)
            attname_mod(dot_pos(j)) = '_';
        end
    end
    
    if ~isempty(white_pos)
        for j = 1:length(white_pos)
            attname_mod(white_pos(j)) = '_';
        end
    end
    
    if ~isempty(point_pos)
        for j = 1:length(point_pos)
            attname_mod(point_pos(j)) = '';
        end
    end
    

    % Write the attribute values to the output structure
    otpt.DataInfo.(attname_mod) = netcdf.getAtt(ncid, globID, attname);
end


% Let's first read all variables to get their names
% THIS IS STILL QUITE NASTY, AS WE HAVE TO CALL netcdf.inqVar and 
% netcdf.inqDim TWICE!!!
if ~strcmp(vars, 'all')
    for i = 1:numvars
        [var_names{i}, xtype(i), dimids{i}, natts] = ...
                                                  netcdf.inqVar(ncid, i-1);
        var_id(i) = i-1;
    end

    req_dims = [];
    for i = 1:length(vars)
        var_indx    = find(strcmp(vars{i}, var_names));
        req_vars(i) = var_id(var_indx);
        req_dims    = unique([req_dims, dimids{var_indx}]);
    end
    
    for i = 1:length(req_dims)
        [dimname{i}, dimlen(i)] = netcdf.inqDim(ncid, req_dims(i));  
        dim_var_indx = strcmp(dimname{i}, var_names);
        req_vars     = unique([req_vars, var_id(dim_var_indx)]);
    end
else
    req_vars = 0:numvars - 1;
    req_dims = 0:numdims - 1;
end


%--------------------------------------------------------------------------
%                               DIMENSIONS
%--------------------------------------------------------------------------

% Loop over the dimensions
for i = 1:length(req_dims)
    % Get the different dimensions of the file
	[dimname{i}, dimlen] = netcdf.inqDim(ncid, req_dims(i));
    % Check for the time-dimension in the file
    if strcmp(dimname{i}, 'time')
        time_id = req_dims(i);
        % Check if a desired time period is provided
        if ~isempty(period)
            dimlen  = period(2);
        end
    elseif strcmp(dimname{i}, 'lon')
        lon_id = req_dims(i);
        % Check if a desired bounding box is provided
        if ~isempty(latlons)
            dimlen = latlons(3);
        end    
    elseif strcmp(dimname{i}, 'lat')
        lat_id = req_dims(i);
        % Check if a desired bounding box is provided
        if ~isempty(latlons)
            dimlen = latlons(4);
        end
    end
        
	% Write the dimensions to the output structure
    if i-1 == unlimdimid
        % Get the unlimited dimension
        otpt.Dimensions.(dimname{i}) = Inf;
    else
        otpt.Dimensions.(dimname{i}) = dimlen;
    end
end


%--------------------------------------------------------------------------
%                                VARIABLES
%--------------------------------------------------------------------------
% Loop over the variables
for i = 1:length(req_vars)
    
    % Get information about the different variables
    [name, xtype, dimids, natts] = netcdf.inqVar(ncid, req_vars(i));

    % Create an empty structure for each variable
    otpt.Variables.(name) = struct();
    
    mval_nc = [];   
    
    % Loop over the attributes
    for j = 1:natts
        % Get the attribute name for the corresponding ID
        att_name = netcdf.inqAttName(ncid, req_vars(i), j-1);   
        
        if strcmp(att_name, '_FillValue')
            % Workaround: As Matlab does not allow a structure-name which
            % starts with an underscore, we simply rename the
            % _FillValue-attribute into FillValue
            mval_nc = netcdf.getAtt(ncid, req_vars(i), att_name);
            if nan_mval == true
                otpt.Variables.(name).FillValue = NaN; 
            else
                otpt.Variables.(name).FillValue = mval_nc;
            end
        elseif ~strcmp(att_name, '_ChunkSizes')
            % Check if the attribute starts with a dash
            if strcmp(att_name(1), '_')
                att_name_new = att_name(2:end);
                
                disp(['Renaming attribute ', att_name, ' to ', ...
                                                             att_name_new])
            else
                att_name_new = att_name;
            end
            
            % Write the other attribute values (except for ChunkSizes; we 
            % don't need that...) to the output structure
            otpt.Variables.(name).(att_name_new) = netcdf.getAtt(ncid, ...
                                                    req_vars(i), att_name);   

        end
    end
    
    % IF the file contains both a _FillValue- and missing_value attribute,
    % it is convenient to remove the missing_value.
    if isfield(otpt.Variables.(name), 'FillValue') && ...
                            isfield(otpt.Variables.(name), 'missing_value')
    	otpt.Variables.(name) = rmfield(otpt.Variables.(name), ...
                                                          'missing_value');
    end

    % Read the data of the variable (but only if the dimids-vector is not
    % empty --> the variable contains some data).
    if ~isempty(dimids)
        % Loop over the variable's dimensions
        hastime = 0;
        haslat  = 0;
        haslon  = 0;
        
        count   = ones(1, length(dimids));
        for j = 1:length(dimids)
            % Get name and length of the dimension
            [dimname, dimlen] = netcdf.inqDim(ncid, dimids(j));
            % Save the name of the dimensions in the dimensions-attribute
            otpt.Variables.(name).dimensions{1, length(dimids)-j+1} = ...
                                                                   dimname;
            % Save the length of each of the variable's dimension  
            count(j) = dimlen;
            if strcmp(dimname, 'time')
                hastime  = 1;
                time_pos = find(dimids == time_id);
            end
            
            if strcmp(dimname, 'lon')
                haslon  = 1;
                lon_pos = find(dimids == lon_id);
            end      
            
            if strcmp(dimname, 'lat')
                haslat = 1;
                lat_pos = find(dimids == lat_id);
            end
            
       
        end      
        
        start = zeros(1, length(dimids));

        if ~isempty(period) & hastime == 1
            start(time_pos) = period(1) - 1;
            count(time_pos) = period(2);
        end
        
        if ~isempty(latlons) & haslon == 1
            start(lon_pos) = latlons(1) - 1;
            count(lon_pos) = latlons(3);
        end
        
        if ~isempty(latlons) & haslat == 1
            start(lat_pos) = latlons(2) - 1;
            count(lat_pos) = latlons(4);
        end
       
        
 
        % Read the data of the actual variable
        tmp = netcdf.getVar(ncid, req_vars(i), start, count);
            
        % Workaround: Switch the dimensions in the data to be conform with
        % the matlab ordering 
        if length(dimids) > 1 && xtype ~= 2
            perm = 1:length(dimids);
            otpt.Data.(name) = permute(tmp, fliplr(perm));
        elseif xtype == 2
            otpt.Data.(name) = cellstr(tmp');
        else
            otpt.Data.(name) = tmp;
        end

        if nan_mval == true
            if ~isempty(mval_nc)
                if ~isnan(mval_nc)
                    % Set the missing values in the output data to NaN
                    otpt.Data.(name)(otpt.Data.(name) == mval_nc) = NaN;
                end
            end
        end
    end        
end

% Close the netCDF-file
netcdf.close(ncid)

% Finally, let's transform the time-vector (if the data contains a
% time-dimension and the time unit is supported)
if find(ismember(fieldnames(otpt.Variables), 'time'), 1)
    if tme_trafo == true
        [otpt.Data.time, otpt.TimeStamp] = ...
                reldate2absdate(otpt.Data.time, otpt.Variables.time.units);
            
        if isfield(otpt.Variables, 'time_bounds')
            tmp_bnds(:, 1:6) = ...
                           reldate2absdate(otpt.Data.time_bounds(:, 1), ...
                                                otpt.Variables.time.units);
                                            
            tmp_bnds(:, 7:12) = ...
                           reldate2absdate(otpt.Data.time_bounds(:, 2), ...
                                                otpt.Variables.time.units);
            otpt.Data.time_bounds = tmp_bnds;
            otpt.Variables.time_bounds.units = 'yyyy-MM-dd HH:mm:ss';   
        end
            
            
        otpt.Variables.time.units = 'yyyy-MM-dd HH:mm:ss';   
        otpt.TimeStamp            = datenum(otpt.Data.time);
             
    end
end


