function [errs, varargout] = batch_validation(obs, mdl, vars, file_prefix, valid_options)
% This function serves as "main function" of the validation module of the
% toolbox. It basically compares the variable "vars" in the two datasets
% obs and mdl. The different performance metrics are defined in a separate
% structure array "valid_options". Please refer to the function
% set_valid_options.m for a description of all different possibilities.
%
% The general workflow of the function is as follows:
% 1. Check if both datasets cover the same time period. If this is not the
%    case, the function truncates both datasets to the same period.
% 2. The function then computes a spatial mask, where all missing values in
%    both datasets are set to NaN. This mask is then used to ensure the
%    same spatial coverage for the computation of the spatial averages.
% 3. Then, the function goes through different levels of performance
%    metrics, depending on the settings in valid_options. For the anomaly
%    metrics, the function removes the (empirically derived) annual cycle
%    from the data. 
%    3.1  Mean 2D long term error fields (full signal)
%    3.2  Mean 2D long term error fields (anomalies)
%    3.3  Mean 2D seasonal error fields (full signal)
%    3.4  Mean 1D long term errors from 2D fields (full signal)
%    3.5  Mean 1D long term errors from 2D fields (anomalies)
%    3.6  Mean 1D seasonal errors from 2D fields (full signal)
%    3.7  Spatial averages (full signal)
%    3.8  Spatial averages (anomalies)
%    3.9  Seasonal spatial averages 
%    3.10 (Optional) Multi-region spatial average.

%
%--------------------------------------------------------------------------
% Input (required):

%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: trunc_TS.m, setmask.m, structerrs.m
%--------------------------------------------------------------------------


% Bring both datasets to the same time-period
if obs.TimeStamp(1, :) ~= mdl.TimeStamp(1, :) | ...
                             obs.TimeStamp(end, :) ~= mdl.TimeStamp(end, :) 

    strt_yr = max([obs.Data.time(1, 1), mdl.Data.time(1, 1)]);
    end_yr  = min([obs.Data.time(end, 1), mdl.Data.time(end, 1)]);

	obs = trunc_TS(obs, strt_yr, end_yr);
    mdl = trunc_TS(mdl, strt_yr, end_yr);
end

% Compute a mask from both datasets, which ensure that the comparison is 
% carried out over the same region
mask = ones(size(obs.Data.(vars), 2), size(obs.Data.(vars), 3));
mask(isnan(obs.Data.(vars)(1, :, :))) = NaN;
mask(isnan(mdl.Data.(vars)(1, :, :))) = NaN;
mdl = setmask(mdl, mask);
obs = setmask(obs, mask);



% -------------------------------------------------------------------------
%               Mean 2D long term error fields (full signal)
% -------------------------------------------------------------------------
disp('Computing Mean 2D long term error fields (full signal)...')
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
disp('Computing Mean 2D long term error fields (anomalies)...')
error_mtrcs_lt_2d_anom           = fieldnames(valid_options.lt_2d_anom);
selected                         = struct2array(valid_options.lt_2d_anom);
error_mtrcs_lt_2d_anom(selected == 0) = [];

if ~isempty(error_mtrcs_lt_2d_anom)
    % Remove the annual cycle from both datasets
    [obs_anom, obs_ann_cycle] = remsc(obs);
    [mdl_anom, mdl_ann_cycle] = remsc(mdl);
    
    errs.lt_2d_anom = structerrs(obs_anom, mdl_anom, vars, ...
                                           error_mtrcs_lt_2d_anom, 'time');
else
    errs.lt_2d_anom = [];
end
disp('Done!')

% -------------------------------------------------------------------------
%                       Mean 2D seasonal error fields 
% -------------------------------------------------------------------------
disp('Computing Mean 2D seasonal error fields...')
error_mtrcs_ssnl_2d                = fieldnames(valid_options.ssnl_2d);
selected                           = struct2array(valid_options.ssnl_2d);
error_mtrcs_ssnl_2d(selected == 0) = [];

if ~isempty(error_mtrcs_ssnl_2d)
    errs.ssnl_2d = seasonal_errors(obs, mdl, vars, error_mtrcs_ssnl_2d);
end
disp('Done!')           
           
% -------------------------------------------------------------------------
%                 Mean 1D long term errors from 2D fields
% -------------------------------------------------------------------------
disp('Computing mean 1D long term errors from 2D fields...')
error_mtrcs_ssnl_1d                = fieldnames(valid_options.ssnl_1d);
selected                           = struct2array(valid_options.ssnl_1d);
error_mtrcs_ssnl_1d(selected == 0) = [];

if ~isempty(error_mtrcs_ssnl_1d)
    for i = 1:length(error_mtrcs_ssnl_1d)
        for j = 1:4
            if isfield(errs.ssnl_2d{j}, error_mtrcs_ssnl_1d{i})
                errs.ssnl_1d.(error_mtrcs_ssnl_1d{i})(1, j) = ...
                    sqrt(nanmean(errs.ssnl_2d{j}.(error_mtrcs_ssnl_1d{i})(:).^2));
            else
                warning('2D Errors not available!')
            end
        end
    end
end
disp('Done!')

% -------------------------------------------------------------------------
%                 Mean 1D seasonal errors from 2D fields
% -------------------------------------------------------------------------
disp('Computing mean 1D seasonal errors from 2D fields...')
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
disp('Computing area averages (full signal)...')
error_mtrcs_spataverage          = fieldnames(valid_options.spataverage);
selected                         = struct2array(valid_options.spataverage);
error_mtrcs_spataverage(selected == 0) = [];

if ~isempty(error_mtrcs_spataverage)
    obs_ts = spatialaverage(obs);
    mdl_ts = spatialaverage(mdl);

    varargout{1}.obs_ts = obs_ts;
    varargout{1}.mdl_ts = mdl_ts;
    
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
disp('Computing area averages (anomalies)...')
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
    

    varargout{1}.obs_ts_anom = obs_ts_anom;
    varargout{1}.mdl_ts_anom = mdl_ts_anom;
    
    errs.spataverage_anom = structerrs(obs_ts_anom, mdl_ts_anom, vars, ...
                                     error_mtrcs_spataverage_anom, 'time');
end
disp('Done!') 

% -------------------------------------------------------------------------
%                        Create the output files
% -------------------------------------------------------------------------

if ~isempty(error_mtrcs_lt_2d)
    % Create a datastructure for the long term error fields (full signal)
    error_struct = spatial_error_output(obs, vars, error_mtrcs_lt_2d, ...
                                                               errs.lt_2d);
    % Add a source to the datastructure
    error_struct.DataInfo.source = [obs.DataInfo.title, ' vs ', ...
                                                       mdl.DataInfo.title];
    error_struct.DataInfo.title = 'Long term mean errors';

    % Create a filename                                              
    fnme_out = [file_prefix, '_2d_lt_', datestr(now, 'yyyy-mm-dd'), '.nc'];
    % Write the output to a netcdf file
    datastruct2netcdf(error_struct, fnme_out);
end
        
if ~isempty(error_mtrcs_lt_2d)    
    % Create a datastructure for the long term error fields (anomalies)                                                               
    anom_error_struct = spatial_error_output(obs, vars, ...
                                  error_mtrcs_lt_2d_anom, errs.lt_2d_anom);
    % Add a source to the datastructure
    anom_error_struct.DataInfo.source = [obs.DataInfo.title, ' vs ', ...
                                                       mdl.DataInfo.title];
    anom_error_struct.DataInfo.title = 'Long term mean errors (anomalies)';                                               
    % Create a filename      
    fnme_out = [file_prefix, '_2d_lt_anom_', datestr(now, 'yyyy-mm-dd'), ...
                                                                    '.nc'];    
    % Write the output to a netcdf file
    datastruct2netcdf(error_struct, fnme_out);
end

if ~isempty(error_mtrcs_ssnl_2d)  
    for i = 1:4
        if i == 1
            ssn = 'MAM';
        elseif i == 2
            ssn = 'JJA';
        elseif i == 3
            ssn = 'SON';
        elseif i == 4
            ssn = 'DJF';
        end
        
        % Create a datastructure for the seasonal error fields 
        ssnl_error_struct{i} = spatial_error_output(obs, vars, ...
                                     error_mtrcs_ssnl_2d, errs.ssnl_2d{i});   
        % Add a title to the datastructure                         
        ssnl_error_struct{i}.DataInfo.source = [obs.DataInfo.title, ...
                                               ' vs ', mdl.DataInfo.title]; 
        ssnl_error_struct{i}.DataInfo.title = ...
                                             ['Seasonal errors for ', ssn];
        % Create a filename      
        fnme_out = [file_prefix, '_2d_ssnl_ssn_', ssn, '_', ...
                                        datestr(now, 'yyyy-mm-dd'), '.nc'];  
        % Write the output to a netcdf file                           
        datastruct2netcdf(ssnl_error_struct{i}, fnme_out);    
    end
end


single_error_output(errs, [file_prefix, '_single_error_metrics_', ...
                                      datestr(now, 'yyyy-mm-dd'), '.txt']);  


                                 					 
                                            
varargout{1} = error_struct;
varargout{2} = anom_error_struct;
varargout{3} = ssnl_error_struct;
