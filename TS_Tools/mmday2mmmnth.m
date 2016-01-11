function out = mmday2mmmnth(inpt)
% The function converts a datastructure with variables in units of
% [mm/day] to [mm/month] 

% First, use the inpt-structure as basis for the output
out = inpt;

% Get all variables
vars = fieldnames(inpt.Variables);

% Get the "fixed" variables
for i = 1:lenght(vars)
    isfixed(i) = isfixedvar(vars{i});
end

% Remove the fixed variables
vars(isfixed == 1) = [];

for i = 1:length(vars)
    % Check if the input variable contains a unit-attribute
    if isfield(inpt.Variables.(vars{i}), 'units')
        % Check if the input variable is in units of mm/month
        if strcmp(inpt.Variables.(vars{i}).units, 'mm/day') | ...
           strcmp(inpt.Variables.(vars{i}).units, 'mm/dy')
            % Compute the number of days for each month
            nrd = eomday(inpt.Data.time(:, 1), inpt.Data.time(:, 2));
            % Divide the data by the number of days
            out.Data.(vars{i}) = bsxfun(@times, out.Data.(vars{i}), nrd);
            % Change the unit of the variable
            out.Variables.(vars{i}).units = 'mm/month';
        end
    end
end
    


new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                      '; MATLAB TS-Tools: mmday2mmmnth.m'];     
out.DataInfo.history = sprintf([new_hist, ' \n', out.DataInfo.history]);