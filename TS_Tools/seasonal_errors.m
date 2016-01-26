function errors = seasonal_errors(obs, mdl, vars, error_mtrcs)
% This function simply extracts the seasonal time-series from obs and mdl
% and then computes the performance metrics (as defined in error_mtrcs) for
% the variable vars
%--------------------------------------------------------------------------
% INPUT:
% - obs, mdl    Input tastructures
% - vars        String which identifies the variable which should be
%               analyzed
% - error_mtrcs Cell array with a list of error metrics, which should be
%               calculated. Refer to set_valid_options.m and structerrs.m
%               for a complete list of performance metrics.
%--------------------------------------------------------------------------
% OUTPUT:
% - errors      1x4 structure array with errors for each season
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: extract_ssnl_ts.m, structerrs.m
%--------------------------------------------------------------------------  

% First, we extract the seasonal average time-series for each dataset
obs_ssnl = extract_ssnl_ts(obs, 'mean');
mdl_ssnl = extract_ssnl_ts(mdl, 'mean');

% Then, we can compute the different error metrics for each season
for i = 1:4
    errors{i} = structerrs(obs_ssnl{i}, mdl_ssnl{i}, vars, error_mtrcs);
end
