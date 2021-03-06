function error_struct = spatial_error_output(inpt_ref, vars, error_mtrcs, errors)
% The function uses an input datastructure (inpt_ref) and corresponding
% errors and creates a datastructure with metadata, ancillary data, ...
% which is required for e.g. writing the errors to netcdf.
%--------------------------------------------------------------------------
% Input (required):
% - inpt_ref        CF-conform data structure which contains the 
% - vars            Analyzed variable
% - error_mtrcs     Cell array with a list of error metrics, which should 
%                   be transformed into a CF datastructure
% - errors          Structure array with the different error data
%    
% Output:   
% - error_struct    Datastructure with errors   
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: create_datastruct.m, full_error_names.m 
%--------------------------------------------------------------------------

% Create a datastructure for the 2D error fields (full signal)
error_struct = create_datastruct(error_mtrcs, '2d_grids', 'double');

% Get the full error names
error_nms    = full_error_names(error_mtrcs);

% Loop over the different error metrics and set some meta-data
for i = 1:length(error_mtrcs)
    error_struct.Variables.(error_mtrcs{i}).long_name = error_nms{i};
    error_struct.Variables.(error_mtrcs{i}).standard_name = error_mtrcs{i};

    if find(ismember(error_mtrcs{i}, {'rmse', 'mae', 'ae'}))
        % Use the unit from the input variable
        err_unit = inpt_ref.Variables.(vars).units;
    elseif find(ismember(error_mtrcs{i}, {'se', 'mse'}))
        % Use the squared unit from the input variable
        err_unit = ['[', obs.Variables.(vars).units, ']^2'];        
    elseif find(ismember(error_mtrcs{i}, {'re', 'are', 'mare', ...
                 'mare_alt', 'sre', 'msre', 'rmsre', 'nrmse', 'cvrmse', ...
                                                        'cvmae', 'pbias'}))
        % Set the error unit to %
        err_unit = '%';
    else
        % In all other cases, set the error unit to no_units
        err_unit = 'no_units';
    end
    
    % Set the error unit
    error_struct.Variables.(error_mtrcs{i}).units = err_unit;
    % Set the identifier for missing values
    error_struct.Variables.(error_mtrcs{i}).FillValue = NaN;
    % Copy the data to the output structure
    error_struct.Data.(error_mtrcs{i}) = errors.(error_mtrcs{i});
end

% Set lat and lon
error_struct.Data.lat       = inpt_ref.Data.lat;
error_struct.Data.lon       = inpt_ref.Data.lon;
error_struct.Dimensions.lat = length(inpt_ref.Data.lat);
error_struct.Dimensions.lon = length(inpt_ref.Data.lon);

% Add some metadata to the output structure
error_struct.DataInfo.references = ...
                     'Spatial Performance metrics from Matlab GeoTS Tools';
error_struct.DataInfo.institution = ...
                'Christof Lorenz, IMK-IFU Garmisch-Partenkirchen, Germany';
            
error_struct.DataInfo.date_created = datestr(now);                       


