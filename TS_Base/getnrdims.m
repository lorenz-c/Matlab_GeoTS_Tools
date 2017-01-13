function [nrdims, varargin] = getnrdims(inpt, vars)

for i = 1:length(vars)
    if isfield(inpt.Variables.(vars{i}), 'dimensions')
        nrdims(i) = length(inpt.Variables.(vars{i}).dimensions);
        for j = 1:length(nrdims(i))
           dimnames{i, j} = inpt.Variables.(vars{i}).dimensions{j};
        end
    else
        nrdims(i) = 0;
    end
end

varargin{1} = dimnames;