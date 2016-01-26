function [istime, vars] = istimevar(inpt, vars)
% Simple function which checks if the variable var in the datastructure
% inpt is a time-series. If this is true, istime is set to 1.
% Otherwise, istime = 0.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - vars        Cell array with a list of variables, which should be
%               checked. If left empty, the function checks all (non-fixed) 
%               variables.
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

    % Loop over the dimensions of each variable
    tmp = find(ismember(dta_dims, 'time'));
    
    if tmp > 0
        istime(i) = 1;
    else
        istime(i) = 0;
    end
end



