function out = mmmnth2mmday(inpt)
% The function converts a datastructure with variables in units of
% [mm/month] to [mm/day] 

% First, use the inpt-structure as basis for the output
out = inpt;

% Get all variables
vars = fieldnames(inpt.Variables);

% Get the "fixed" variables
for i = 1:length(vars)
    isfixed(i) = isfixedvar(vars{i});
end

% Remove the fixed variables
vars(isfixed == 1) = [];

for i = 1:length(vars)
    % Check if the input variable contains a unit-attribute
    if isfield(inpt.Variables.(vars{i}), 'units')
        % Check if the input variable is in units of mm/month
        if strcmp(inpt.Variables.(vars{i}).units, 'mm/month') | ...
           strcmp(inpt.Variables.(vars{i}).units, 'mm/mnth')
            disp(['mmmnth2mmday.m: Found matching variable: ', vars{i}])
            % Compute the number of days for each month
            nrd = eomday(inpt.Data.time(:, 1), inpt.Data.time(:, 2));
            % Get the "position" of the time-dimension
            dimpos = getdimpos(inpt, vars{i}, 'time');
             % Divide the data by the number of days
            if dimpos == 1
                out.Data.(vars{i}) = bsxfun(@rdivide, ...
                                                out.Data.(vars{i}), nrd);
            elseif dimpos == 2
                out.Data.(vars{i}) = bsxfun(@times, ...
                                                out.Data.(vars{i}), nrd');
            end  
           
            % Change the unit of the variable
            out.Variables.(vars{i}).units = 'mm/day';
        end
    end
end
    


new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                      '; MATLAB TS-Tools: mmmnth2mmday.m'];     
out.DataInfo.history = sprintf([new_hist, ' \n', out.DataInfo.history]);