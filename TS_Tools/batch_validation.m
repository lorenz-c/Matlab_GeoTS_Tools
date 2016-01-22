function errs = batch_validation(obs, mdl, vars, valid_options)
    
% Bring both datasets to the same time-period
if obs.TimeStamp(1, :) ~= mdl.TimeStamp(1, :) | ...
                             obs.TimeStamp(end, :) ~= mdl.TimeStamp(end, :) 

    strt_yr = max([obs.Data.time(1, 1), mdl.Data.time(1, 1)]);
    end_yr  = min([obs.Data.time(end, 1), mdl.Data.time(end, 1)]);

	obs = trunc_TS(obs, strt_yr, end_yr);
    mdl = trunc_TS(mdl, strt_yr, end_yr);
end

% Compute the seasonal averages
obs_ssnl_lt = ts_average(obs, 'seasonal_lt');
obs_mnth_lt = ts_average(obs, 'monthly_lt');

mdl_ssnl_lt = ts_average(mdl, 'seasonal_lt');
mdl_mnth_lt = ts_average(mdl, 'monthly_lt');


% -------------------------------------------------------------------------
%               Mean 2D long term error fields (full signal)
% -------------------------------------------------------------------------
disp('Computing Mean 2D long term error fields (full signal)')
% Get the selected metrics
error_mtrcs_lt_2d                = fieldnames(valid_options.lt_2d);
selected                         = struct2array(valid_options.lt_2d);
error_mtrcs_lt_2d(selected == 0) = [];

if ~isempty(error_mtrcs_lt_2d)
    errs.lt_2d = structerrs(obs, mdl, vars, error_mtrcs_lt_2d, 'time');
else
    errs.lt_2d = [];
end
disp('Done!')

% -------------------------------------------------------------------------
%                Mean 2D long term error fields (anomalies)
% -------------------------------------------------------------------------
disp('Mean 2D long term error fields (anomalies)')
% Get the selected metrics
error_mtrcs_lt_2d_anom           = fieldnames(valid_options.lt_2d_anom);
selected                         = struct2array(valid_options.lt_2d_anom);
error_mtrcs_lt_2d_anom(selected == 0) = [];

if ~isempty(error_mtrcs_lt_2d_anom)
    % Compute the seasonal averages
    obs_anom = remsc(obs);
    mdl_anom = remsc(mdl);
    
    errs.lt_2d_anom = structerrs(obs_anom, mdl_anom, vars, ...
                                           error_mtrcs_lt_2d_anom, 'time');
else
    errs.lt_2d_anom = [];
end
disp('Done!')

% -------------------------------------------------------------------------
%                 Mean 1D long term errors from 2D fields
% -------------------------------------------------------------------------
disp('Computing mean 1D long term errors from 2D fields')
% Get the selected metrics
error_mtrcs_lt_1d                = fieldnames(valid_options.lt_1d);
selected                         = struct2array(valid_options.lt_1d);
error_mtrcs_lt_1d(selected == 0) = [];

if ~isempty(error_mtrcs_lt_1d)
    for i = 1:length(error_mtrcs_lt_1d)
        if isfield(errs.lt_2d, error_mtrcs_lt_1d{i})
            errs.lt_1d.(error_mtrcs_lt_1d{i}) = ...
                    sqrt(nanmean(errs.lt_2d.(error_mtrcs_lt_1d{i})(:).^2));
        else
            tmp = structerrs(obs, mdl, vars, error_mtrcs_lt_1d{i}, 'time');
            errs.lt_1d.(error_mtrcs_lt_1d{i}) = ...
                        sqrt(nanmean(tmp.(error_mtrcs_lt_1d{i})(:).^2));
        end
    end
end
disp('Done!')
    
% -------------------------------------------------------------------------
%                       Area averages (full signal)
% -------------------------------------------------------------------------
disp('Computing area averages (full signal)')
error_mtrcs_spataverage          = fieldnames(valid_options.spataverage);
selected                         = struct2array(valid_options.spataverage);
error_mtrcs_spataverage(selected == 0) = [];

if ~isempty(error_mtrcs_spataverage)
    obs_ts = spatialaverage(obs);
    mdl_ts = spatialaverage(mdl);
    
    errs.spataverage = structerrs(obs_ts, mdl_ts, vars, ...
                                          error_mtrcs_spataverage, 'time');
else
    obs_ts = [];
    mdl_ts = [];
end
disp('Done!')   
    
% -------------------------------------------------------------------------
%                       Area averages (anomalies)
% -------------------------------------------------------------------------
disp('Computing area averages (anomalies)')
error_mtrcs_spataverage_anom     = ...
                                fieldnames(valid_options.spataverage_anom);
selected                         = ...
                              struct2array(valid_options.spataverage_anom);
error_mtrcs_spataverage_anom(selected == 0) = [];

if ~isempty(error_mtrcs_spataverage_anom)
    if isempty(obs_ts) 
        obs_ts = spatialaverage(obs);
    end
    if isempty(mdl_ts)
        mdl_ts = spatialaverage(mdl);
    end
    
    obs_ts_anom = remsc(obs_ts);
    mdl_ts_anom = remsc(mdl_ts);
    
    errs.spataverage_anom = structerrs(obs_ts_anom, mdl_ts_anom, vars, ...
                                     error_mtrcs_spataverage_anom, 'time');
end
disp('Done!') 

