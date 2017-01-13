function [istime, vars] = istimevar(inpt, vars)
% Simple function which checks if the variable var in the datastructure
% inpt is a time-series. If this is true, istime is set to 1.
% Otherwise, istime = 0.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - vars        Cell array with a list of variables, which should be
%               checked. If left empty, the function checks all variables.
%--------------------------------------------------------------------------
% OUTPUT:
% - istime      Logical variable: 1 -> variable has time dimension
%                                 0 -> variable has no time dimension
% - vars        List of variables, which have been checked
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------   

if nargin < 2, vars = fieldnames(inpt.Variables); end

% Get the data dimensions
for i = 1:length(vars)
    if isfield(inpt.Variables.(vars{i}), 'dimensions')
        dta_dims = inpt.Variables.(vars{i}).dimensions;

        % Search for "time" in the list of dimensions
        istime(i) = ismember('time', dta_dims);
    else
        istime(i) = 0;
    end
end




