function [isgrid, vars] = isgridvar(inpt, vars)
% Simple function which checks if the variable var in the datastructure
% inpt has lat- and lon-coordinates. If this is true, isgrid is set to 1.
% Otherwise, isgrid = 0.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - vars        Cell array with a list of variables, which should be
%               checked. If left empty, the function checks all variables.
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
% Uses: 
%--------------------------------------------------------------------------   

if nargin < 2, vars = fieldnames(inpt.Variables); end


for i = 1:length(vars)
    % Get the data dimensions
    if isfield(inpt.Variables.(vars{i}), 'dimensions')
        dta_dims = inpt.Variables.(vars{i}).dimensions;
  
        has_lat = max(ismember({'lat', 'latitude'}, dta_dims));
        has_lon = max(ismember({'lon', 'longitude'}, dta_dims));
    
        has_x   = ismember({'x'}, dta_dims);
        has_y   = ismember({'y'}, dta_dims);
        
        if has_lat == 1 && has_lon == 1
            isgrid(i) = 1;
        elseif has_x == 1 && has_y == 1
            isgrid(i) = 2;
        else
            isgrid(i) = 0;
        end
    else
        isgrid(i) = 0;
    end
end




