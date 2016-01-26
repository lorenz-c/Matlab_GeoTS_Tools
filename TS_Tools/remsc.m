function [Res, Ssnl_cycle] = remsc_ts(inpt, vars)




if nargin < 2, vars = 'all'; end

if strcmp(vars, 'all')
    vars = fieldnames(inpt.Variables);
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
    end
    
    % Remove all fixed variables
    vars(isfixed == 1) = [];
end



% The output will have the same "layout" like the input
Res = inpt;

% Compute the seasonal cycle
Ssnl_cycle = ts_average(inpt, 'monthly_lt');

for i = 1:length(vars)
    % Get the position of the time dimension
    dta_dims = inpt.Variables.(vars{i}).dimensions;
    dimpos   = find(ismember(dta_dims, 'time'));
    
    for j = 1:12
        % Get the indices for each month
        mnth_indx = find(inpt.Data.time(:, 2) == j);

        if length(dta_dims) <= 2
            if dimpos == 1
                Res.Data.(vars{i})(mnth_indx, :) = bsxfun(@minus, ...
                                      Res.Data.(vars{i})(mnth_indx, :), ...
                                          Ssnl_cycle.Data.(vars{i})(j, :));
            elseif dimpos == 2
                Res.Data.(vars{i})(:, mnth_indx) = bsxfun(@minus, ...
                                      Res.Data.(vars{i})(:, mnth_indx), ...
                                          Ssnl_cycle.Data.(vars{i})(:, j));
            end
        elseif length(dta_dims) == 3
            if dimpos == 1
                Res.Data.(vars{i})(mnth_indx, :, :) = bsxfun(@minus, ...
                                   Res.Data.(vars{i})(mnth_indx, :, :), ...
                                       Ssnl_cycle.Data.(vars{i})(j, :, :));
            else
                error('First dimension must be "time" for 3D arrays!')
            end

        end
    end
 
    % Add the applied method to the source-attribute
    new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                              ': Seasonal cycle removed '];  
    if isfield(Res.Data.(vars{i}), 'source')
        Res.Variables.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                          Res.Variables.(vars{i}).source]);
    else
        Res.Variables.(vars{i}).source = new_srce;
    end
    
    
end
