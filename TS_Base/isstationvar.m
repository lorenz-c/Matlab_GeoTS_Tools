function [isstation, vars] = isstationvar(inpt, vars)
% Simple function which checks if the variable var in the datastructure
% inpt is station data. If this is true, isstation is set to 1.
% Otherwise, isstation = 0.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - vars        Cell array with a list of variables, which should be
%               checked. If left empty, the function checks all variables.
%--------------------------------------------------------------------------
% OUTPUT:
% - isstation   Logical variable: 1 -> variable has station dimension
%                                 0 -> variable has no station dimension
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
         
for i = 1:length(vars)
    % Get the data dimensions
    dta_dims = inpt.Variables.(vars{i}).dimensions;
    % Search for "time" in the list of dimensions
    isstation(i) = ismember('stations', dta_dims);
end



