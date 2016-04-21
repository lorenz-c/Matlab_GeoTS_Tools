function otpt = setmask(inpt, bbox)

% Copy the input to the output
otpt = inpt;

% Get the variables of the input data
vars = fieldnames(inpt.Variables);
for i = 1:length(vars)
    isfixed(i) = isfixedvar(vars{i});
end
% Remove the fixed variables
vars(isfixed == 1) = [];

% Check for gridded variables
isgrid = isgridvar(inpt, vars);
% Remove variables which do not have the lat- and lon-dimension
vars(isgrid == 0) = [];

for i = 1:length(vars)
    % Get the dimension
    dims = inpt.Variables.(vars{i}).dimensions;
    
    if length(dims) == 3
       % keyboard
        mask_3d                = reshape(mask, [1 size(mask)]);

        otpt.Data.(vars{i}) = bsxfun(@times, otpt.Data.(vars{i}), mask_3d);
    elseif length(dims) == 2
        otpt.Data.(vars{i}) = otpt.Data.(vars{i}).*mask;
    end
end

% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                           '; MATLAB TS-Tools: setmask.m'];
if isfield(inpt.DataInfo, 'history')           
    otpt.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                   inpt.DataInfo.history]);
else
    otpt.DataInfo.history = new_hist;
end

