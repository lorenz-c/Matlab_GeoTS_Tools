function ssnl_ts = extract_ssnl_ts(inpt, method, vars)
% This function extracts four time-series (one for each season) from an
% input datastructure inpt. Therefore, it can be chosen if all time-steps,
% which correspont to a single season, are extracted (method = 'full') or
% if some averaging should be performed. Please have a look at the
% ts_average.m function for all different options.
% If no variables are selected, the function automatically uses all
% variables with a time dimension and extracts the seasonal time series.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Input datastructure
% - method      Can be set to 'full' (default) -> All time-steps which
%               correspond to a single season are extracted. Other options
%               allow the averaging of the seasons. Therefore, method can
%               be set to (e.g.) 'mean', 'sum', ...
%               Refer to the documentation of ts_average.m for further
%               options.
% - vars        By default, the function extracts the seasonal time-series
%               for all variables with a time dimension.
%--------------------------------------------------------------------------
% OUTPUT:
% - ssnl_ts     1x4 structure array, which contains the seasonal
%               time-series for each of the four seasons in the following
%               order: 
%               ssnl_ts{1} -> MAM 
%               ssnl_ts{2} -> JJA 
%               ssnl_ts{3} -> SON
%               ssnl_ts{4} -> DJF
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: isfixedvar.m, istimevar.m, ts_average.m, getdimpos.m
%--------------------------------------------------------------------------   
if nargin < 3, vars = 'all'; end
if nargin < 2, method = 'full'; end


if strcmp(vars, 'all')
    % Remove fixed variables
    vars               = fieldnames(inpt.Variables);
    isfixed            = isfixedvar(vars);
    vars(isfixed == 1) = [];
    
    % Check for variables with a "time"-dimension
    istime            = istimevar(inpt, vars);
    vars(istime == 0) = [];
end

if strcmp(method, 'full')
    % Get a vector with all months of the data
    mnths = inpt.Data.time(:, 2);
    % Go through the four seasons and save the indices of the corresponding
    % time-steps
    for i = 1:4
        if i == 4
            indx_1 = find(mnths == 12);
            indx_2 = find(mnths == 1);
            indx_3 = find(mnths == 2);            
        else
            indx_1 = find(mnths == (i+1)*3-3);
            indx_2 = find(mnths == (i+1)*3-2);
            indx_3 = find(mnths == (i+1)*3-1);
        end
        % Save all indices from a single season in tmp
        tmp     = [indx_1; indx_2; indx_3];
        % Order tmp (ascending) and save the result in indx
        indx{i} = sort(tmp, 'ascend');
    end
       
else
    % Compute the seasonal average (or sum, ...)
    inpt  = ts_average(inpt, 'seasonal', method);
    % Get all months
    mnths = inpt.Data.time(:, 2);   
    % Get the indices of the different seasons. Note that (as we've already
    % computed an average) the seasons are now identified by the months 4,
    % 7, 10, and 1, respectively!
    indx{1} = find(mnths == 4);
    indx{2} = find(mnths == 7);
    indx{3} = find(mnths == 10);
    indx{4} = find(mnths == 1);
end


            
for i = 1:4
    % Use inpt as a "template"
    ssnl_ts{i} = inpt;
    
    % Get the coorect time-data
    ssnl_ts{i}.Data.time = inpt.Data.time(indx{i}, :);
    ssnl_ts{i}.TimeStamp = inpt.TimeStamp(indx{i});
    
    % Loop over the variables
    for j = 1:length(vars)
        % Get the position of the time dimension for each variable
        dimpos = getdimpos(inpt, vars{j}, 'time');
        % Get the number of dimensions for each variable
        dimlnth = length(inpt.Variables.(vars{j}).dimensions); 
        % Depending on the lenght of the variable's dimension array and the
        % "position" of the time dimension, we now extract all data
        % according to the indices in indx{i}.
    	if dimlnth <= 2 
            if dimpos == 1
                ssnl_ts{i}.Data.(vars{j}) = ...
                                           inpt.Data.(vars{j})(indx{i}, :);
            elseif dimpos == 2
                ssnl_ts{i}.Data.(vars{j}) = ...
                                           inpt.Data.(vars{j})(:, indx{i});
            end
        elseif dimlnth == 3
            if dimpos == 1
                ssnl_ts{i}.Data.(vars{j}) = ...
                                        inpt.Data.(vars{j})(indx{i}, :, :);
            elseif dimpos == 2
                ssnl_ts{i}.Data.(vars{j}) = ...
                                        inpt.Data.(vars{j})(:, indx{i}, :);
            elseif dimpos == 3
                ssnl_ts{i}.Data.(vars{j}) = ...
                                        inpt.Data.(vars{j})(:, :, indx{i});
            end
        else
            error('Unknown shape of the data')
        end

        % Add some short description of the calculations to the source-
        % attribute of each variable
        new_srce = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                         ': Data for season ', num2str(i)];   
        if isfield(ssnl_ts{i}.Variables.(vars{j}), 'source')
            ssnl_ts{i}.Variables.(vars{j}).source = sprintf([new_srce, ...
                            ' \n', ssnl_ts{i}.Variables.(vars{j}).source]);
        else
            ssnl_ts{i}.Variables.(vars{j}).source = new_srce;
        end   
    end
    
    % Update the file history
    new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                   '; MATLAB TS-Tools: extract_ssnl_ts.m'];    
        
    ssnl_ts{i}.DataInfo.history = ...
                   sprintf([new_hist, ' \n', ssnl_ts{i}.DataInfo.history]);
end       

        
    
   
    
 

            