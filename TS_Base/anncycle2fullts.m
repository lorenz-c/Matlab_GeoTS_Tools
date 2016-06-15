function otpt = anncycle2fullts(ann_cycle, dtes, vars)
% The function re-sizes the (climatological) annual cycle (which has only
% 12 time steps) to the length of a given time vector dtes. 
%--------------------------------------------------------------------------
% INPUT:
% - ann_cycle   Datastructure, which contains the long-term annual cycle
%               (e.g. from the ts_average.m function)
% - dtes        Nx6 matrix, which contains yyyy-mm-dd hh:mm:ss
% - vars        Cell array with a list of variables, which should be
%               re-scaled
%--------------------------------------------------------------------------
% OUTPUT:
% - otpt        Datastructure where each variable in vars has now the
%               length of dtes.
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: isfixedvar.m, istimevar.m
%--------------------------------------------------------------------------

if nargin < 3, vars = 'all'; end

otpt           = ann_cycle;
otpt.Data.time = dtes;
otpt.TimeStamp = datenum(dtes);


if strcmp(vars, 'all')
    % First, remove the fixed variables
    vars               = fieldnames(ann_cycle.Variables);
    isfixed            = isfixedvar(vars);
    vars(isfixed == 1) = [];
    
    % Then, remove the variables without time dimension
    istime             = istimevar(ann_cycle, vars);
    vars(istime == 0)  = [];
end


for i = 1:length(vars)
    % Check if the current variable is found in the datastructure
    if ismember(vars{i}, fieldnames(ann_cycle.Variables))
        % Get all dimensions of the current variable
        dta_dims = ann_cycle.Variables.(vars{i}).dimensions;
        % Get the "position" of the time dimension
        dimpos   = find(ismember(dta_dims, 'time'));
        % Loop over twelve months
        for j = 1:12
            % Get the indices for each month
            mnth_indx = find(dtes(:, 2) == j);
        
            if length(dta_dims) <= 2
                % Variable contains one (or more) discrete time-series
                if dimpos == 1 
                    % Time -> First dimension
                    otpt.Data.(vars{i})(mnth_indx, :) = ...
                                            ann_cycle.Data.(vars{i})(j, :);
                elseif dimpos == 2
                    % Time -> Second dimension
                    otpt.Data.(vars{i})(:, mnth_indx) = ...
                                            ann_cycle.Data.(vars{i})(:, j);
                end
            elseif length(dta_dims) == 3
                % Variable contains gridded time-series
                if dimpos == 1
                    % Time -> First dimension
                    otpt.Data.(vars{i})(mnth_indx, :, :) = ...
                              repmat(ann_cycle.Data.(vars{i})(j, :, :), ...
                                                  length(mnth_indx), 1, 1);
                else
                    error('First dimension must be "time" for 3D arrays!')
                end                                                                     
            elseif length(dta_dims) == 4
                % Variable contains gridded time-series
                if dimpos == 1
                    % Time -> First dimension
                    otpt.Data.(vars{i})(mnth_indx, :, :, :) = ...
                           repmat(ann_cycle.Data.(vars{i})(j, :, :, :), ...
                                               length(mnth_indx), 1, 1, 1);
                else
                    error('First dimension must be "time" for 3D arrays!')
                end
            end
        end
    else
        warning(['Variable ', vars{i}, ' not found in datastructure'])
    end
end

% Update the file history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                  '; MATLAB TS-Tools: anncycletofullts.m'];
        
otpt.DataInfo.history = sprintf([new_hist, ' \n', otpt.DataInfo.history]);        
    
    
    
    


