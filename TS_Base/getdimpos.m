function dimpos = getdimpos(inpt, varnames, dimnames)
% This small function simply gives the indices (i.e. the "positions") of a
% list of dimensions "dimnames" for the variables "varnames" in the 
% datastructure "inpt".
%--------------------------------------------------------------------------
% Input (required):
% - inpt        CF-conform data structure
% - varnames    Cell array or string with variable names
% - dimnames    Cell array or string with dimension names

% Output:
% - dimpos      Matrix, vector, or scalar (depending on the input
%               parameter) with the dimension index for each variable and 
%               dimension
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Transform the string in a cell array if only a single variable is 
% provided.
if isstr(varnames)
    vars{1} = varnames;
else
    vars    = varnames;
end

% Transform the string in a cell array if only a single dimension is 
% provided.
if isstr(dimnames)
    dims{1} = dimnames;
else
    dims    = dimnames;
end

% Loop over the variables
for i = 1:length(vars)
    % Get the dimensions of the variable
    dims_var  = inpt.Variables.(vars{i}).dimensions;
    
    % Loop over the dimensions
    for j = 1:length(dims)
        % Check if dims{j} appears in the variable's dimensions
        ismem = ismember(dims_var, dims{j});
    
        % Get the position of the dimension
        tmp = find(ismem == 1);
        
        if isempty(tmp)
            % If the actual dimension is not found in the variable, the
            % function gives a warning and the index is set to NaN.
            warning(['Dimension ', dims{j}, ' not found in ', vars{i}]);
            dimpos(i, j) = NaN;
        else
            % Else, the index is added to the ouput
            dimpos(i, j) = find(ismem == 1);
        end
    end
end




