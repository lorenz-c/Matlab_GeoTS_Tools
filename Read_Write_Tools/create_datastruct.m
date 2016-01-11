function out = create_datastruct(varnames, vartypes, varprec)
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
%               '3d_grids', and '4d_grids'. TBA: Profiles
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

% 1. Let's create an empty datastructure with all mandatory fields
out = struct('DataInfo', [], ...
             'Dimensions', [], ...
             'Variables', [], ...
             'Data', []);
                                                           
% 2. In this function, we only create the mandatory global metadata
out.DataInfo = struct('title', char.empty(0), ...
    'source', char.empty(0), ...
    'history', [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                           '; create_CF_datastruct.m'], ...
    'references', char.empty(0), ...
    'comment', char.empty(0), ...
    'Conventions', 'CF-1.6', ...
    'institution', char.empty(0), ...
    'cdm_data_type', char.empty(0));

% 3. Add the variables by calling the addvariable.m-function
out = addvariable(out, varnames, vartypes, varprec)