function isvar = isfixedvar(varnme)
% The function checks if the variable given by varname is a fixed variable.
% Currently, there are the following fixed variables:
% 'latitude', 'longitude', 'lat', 'lon', 'time', 'Time', 'times', 'Times',
% 'regions', 'region_map', 'region_names', 'levels', 'station_name',
% 'stations', 'alt'.
%--------------------------------------------------------------------------
% INPUT:
% - varnme      Variable name (string) which should be checked
%--------------------------------------------------------------------------
% OUTPUT:
% - isvar       Logical variable: 1 -> variable is a fixed variable
%                                 0 -> variable is no fixed variable
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% List of "fixed" variable names
fixed_vars = {'latitude', ...
              'longitude', ...
              'lat', ...
              'lon', ...
              'time', ...
              'Time', ...
              'times', ...
              'Times', ...
              'regions', ...
              'region_map', ...
              'region_names', ...
              'region_area', ...
              'levels', ...
              'station_name', ...
              'stations', ...
              'alt', ...
              'climatology_bounds', ...
              'time_bounds'};
          
% Check each variable in varnme
isvar = ismember(varnme, fixed_vars);

          


