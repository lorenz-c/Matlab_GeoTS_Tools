function [vars, dimvars] = remdims(inpt)
% The function extracts a list of variables of a datastructure inpt, which
% are not used as dimension variables (such as, e.g., time, lat, lon).
%--------------------------------------------------------------------------
% Input (required):
% - inpt        Datastructure with one or more variables, which are not
%               dimension variables
% Output:
% - vars        List of variables, which are not used as dimension
%               variables
% - dimvars     List of dimension variables
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Dimensions of the data
dims = fieldnames(inpt.Dimensions);

% Get the names of the variables
vars = fieldnames(inpt.Variables);

% Check if the dimensions also have a corresponding variable
for i = 1:length(dims)
    tmp = find(ismember(vars, dims{i}));
    if isempty(tmp)
        dim_ids(i) = NaN;
    else
        dim_ids(i) = tmp;
    end
end

% Remove the dimensions which do not have a corresponding variable
dim_ids(isnan(dim_ids)) = [];

% Extract the dimension variables
dimvars = vars(dim_ids, :);

% Remove the dimension variables
vars(dim_ids, :) = [];