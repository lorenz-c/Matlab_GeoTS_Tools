function [isgrid, vars] = isgridvar(inpt, vars)
% Simple function which checks if the variable var in the datastructure
% inpt has lat- and lon-coordinates. If this is true, isgrid is set to 1.
% Otherwise, isgrid = 0.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - vars        Cell array with a list of variables, which should be
%               checked. If left empty, the function checks all (non-fixed) 
%               variables.
%--------------------------------------------------------------------------
% OUTPUT:
% - isgrid      Logical variable: 1 -> variable is a gridded variable
%                                 0 -> variable is no gridded variable
% - vars        List of variables, which have been checked
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: isfixedvar.m
%--------------------------------------------------------------------------   

if nargin < 2, vars = 'all'; end

if strcmp(vars, 'all')
    vars = fieldnames(inpt.Variables);
    
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
    end
    
    % Remove the fixed variables from the data
    vars(isfixed == 1) = [];
end
         
for i = 1:length(vars)
    % Get the data dimensions
    dta_dims = inpt.Variables.(vars{i}).dimensions;

    has_lat = 0;
    has_lon = 0;
    
    % Loop over the dimensions of each variable
    for j = 1:length(dta_dims)
        % Check for latitude
        if strcmp(dta_dims{j}, 'lat') | strcmp(dta_dims{j}, 'latitude')
            has_lat = 1;
        end
        % Check for longitude
        if strcmp(dta_dims{j}, 'lon') | strcmp(dta_dims{j}, 'longitude')
            has_lon = 1;
        end
    end

	if has_lat == 1 && has_lon == 1
        isgrid(i) = true;
    else
        isgrid(i) = false;
    end
end




