function [otpt] = netcdf2datastruct(fnme, tme_trafo, nan_mval, vars)
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
%               TBA: Support for time zones (e.g, UTC +01)
% - nan_mval    Boolean variable: if set to true, the function converts the
%               missing values into NaNs.
%               provided, the _FillValue in the netCDF-file is used as 
%               identifier. However, it is more convenient to work with 
%               matlab's NaNs.
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
    dash_pos = find(ismember(attname, '-'));
    dot_pos  = find(ismember(attname, '.'));
    
    if ~isempty(dash_pos)
        attname_new = attname;
        for j = 1:length(dash_pos)
            attname_new(dash_pos(j)) = '_';
        end
    else
        attname_new = attname;
    end
    
    if ~isempty(dot_pos)
        attname_new = attname;
        for j = 1:length(dot_pos)
            attname_new(dot_pos(j)) = '_';
        end
    else
        attname_new = attname;
    end
    
    % Write the attribute values to the output structure
    otpt.DataInfo.(attname_new) = netcdf.getAtt(ncid, globID, attname);
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
        elseif strcmp(att_name, '_ChunkSizes')
            otpt.Variables.(name).ChunkSizes = netcdf.getAtt(ncid, ...
                                                    req_vars(i), att_name);    
        else
            % Write the attribute values to the output structure
            otpt.Variables.(name).(att_name) = netcdf.getAtt(ncid, ...
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
        for j = 1:length(dimids)
            % Get name and length of the dimension
            [dimname, dimlen] = netcdf.inqDim(ncid, dimids(j));
            % Save the name of the dimensions in the dimensions-attribute
            otpt.Variables.(name).dimensions{1, length(dimids)-j+1} = ...
                                                                   dimname;
        end      
        
        % Read the data of the actual variable
        tmp = netcdf.getVar(ncid, req_vars(i));

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
if isfield(otpt, 'DataInfo')
    if isfield(otpt.DataInfo, 'geospatial_lat_resolution')
        otpt.DataInfo.geospatial_lat_resolution = abs(otpt.Data.lat(2) - ...
                                                         otpt.Data.lat(1));
    end
    if isfield(otpt.DataInfo, 'geospatial_lon_resolution')
        otpt.DataInfo.geospatial_lon_resolution = abs(otpt.Data.lon(2) - ...
                                                         otpt.Data.lon(1));
    end
end


% Finally, let's transform the time-vector (if the data contains a
% time-dimension)
if find(ismember(fieldnames(otpt.Variables), 'time'), 1)
    if tme_trafo == true
        [otpt.Data.time, otpt.TimeStamp] = ...
                reldate2absdate(otpt.Data.time, otpt.Variables.time.units);
        otpt.Variables.time.units = 'yyyy-MM-dd HH:mm:ss';   
        otpt.TimeStamp            = datenum(otpt.Data.time);
    end
end




