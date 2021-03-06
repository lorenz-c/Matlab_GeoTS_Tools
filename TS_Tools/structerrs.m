function [errs, varargout] = structerrs(obs_struct, mdl_struct, vars, err_mtrcs, dim, params, val_period, anom)
%--------------------------------------------------------------------------
% The function computes error metrics between two datastructures. It allows
% to select different error quantities and the dimension along which the
% errors should be calculated.
%
% Currently, the following metrics are supported:
% - ae        Absolute errors 
% - mae       Mean absolute errors
% - se        Squared errors
% - mse       Mean squared errors
% - rmse      Root mean squared errors
% - re        Relative errors (w.r.t. the data from struct1 --> obs)
% - mre       Absolute relative errors (abs(re))
% - mare      Mean absolute relative errors
% - mare_2    Mean absolute relative errors (w.r.t. the maximal amplitude 
%             of the data from struct1)
% - sre       Squared relative errors
% - msre      Mean squared relative errors
% - rmsre     Root mean squared relative errors
% - nrmsre    Normalized (w.r.t. the maximal amplitude of the data from 
%             struct1) root mean squared relative error
% - cvrmse    Coefficient of variation of the RMSE
% - cvmae     Coefficient of variation of the MAE
% - cod       Coefficient of determination
% - pbias     Relative bias (w.r.t. struct1)
% - corr      Pearsons correlation coefficient
% - cov_0     Covariance (normalized by N)
% - cov_1     Covariance (normalized by N - 1)
% 
%
% Depending on the chosen error metric, the output will have the same
% dimensions as the input data. If 'average' error metrics are required,
% the dimension along which the errors are computed is removed.
%
% The function already contains some placeholders for more metrics, which
% will be successively added.
%--------------------------------------------------------------------------
% INPUT:
% - obs       CF datastructure (e.g. with observations)
% - mdl       CF datastructure
% - vars      Variable for which the errors should be computed 
%             (e.g. 'prec')
% - err_mtrcs Cell-array of error metrics, which should be calculated
% - dim       Dimension along which the errors should be computed (e.g.
%             'time', 'lat', 'stations', ...)
% - params    Parameters for some error metrics (KGE, ...)
%--------------------------------------------------------------------------
% OUTPUT:
% - errs      Vector (or array) with errors. 
% - varargout The first element contains the time period, over which the
%             errors are computed.
%             The second element has the same size as errs and contains the
%             number of valid datapoints, which were used for the
%             calculation.
%--------------------------------------------------------------------------
% EXAMPLE:
% >> [errs, tme_out, nr_dta] = structerrs(dta1, dta2, 'prec', 'mae', ...
%                               'time');
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January 2016
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: trunc_TS.m
%--------------------------------------------------------------------------
% References: 
% - P. Krause, D. P. Boyle, and F. B?se, 2006: Comparison of different
%   efficiency criteria for hydrological model assessment, Advances in
%   Geosciences, 5, 89 - 97
%--------------------------------------------------------------------------

if nargin < 8, anom       = false;    end
if nargin < 7, val_period = -1;       end
if nargin < 6, params     = [1 1 1];  end
if nargin < 5, dim        = 'time';   end
if nargin < 4, err_mtrcs  = {'rmse'}; end

% First, check if the variables from both datasets have the same unit
% if isfield(obs_struct.Variables.(vars{1}), 'unit') & ...
%                                isfield(mdl_struct.Variables.(vars{2}), 'unit')
%     if obs_struct.Variables.(vars{1}).unit ~= ...
%                                            mdl_struct.Variables.(vars{2}).unit
%         error('Both variables must have the same unit!')
%     end
% end

% Truncate the two structures to the same time period
if val_period ~= -1
    obs_struct = trunc_TS(obs_struct, val_period(1, :), val_period(2, :));
    mdl_struct = trunc_TS(mdl_struct, val_period(1, :), val_period(2, :));
else
    if obs_struct.TimeStamp(1, :) ~= mdl_struct.TimeStamp(1, :) | ...
               obs_struct.TimeStamp(end, :) ~= mdl_struct.TimeStamp(end, :) 

        strt_yr = max([obs_struct.Data.time(1, 1), ...
                                              mdl_struct.Data.time(1, 1)]);
        end_yr  = min([obs_struct.Data.time(end, 1), ...
                                            mdl_struct.Data.time(end, 1)]);

        obs_struct = trunc_TS(obs_struct, strt_yr, end_yr);
        mdl_struct = trunc_TS(mdl_struct, strt_yr, end_yr);
    end
end

if anom == true
    obs_struct = remsc(obs_struct);
    mdl_struct = remsc(mdl_struct);
end

% Save the reference period in tme_out
tme_out.Data      = obs_struct.Data.time;
tme_out.TimeStamp = obs_struct.TimeStamp;

% Get the data for the desired variable from both datasets
obs = obs_struct.Data.(vars{1});
mdl = mdl_struct.Data.(vars{2});

% Get the dimensions of the data
dta_dims = obs_struct.Variables.(vars{1}).dimensions;

% Get the "position" of the dimension, over which the errors should be
% calculated
dimpos = find(ismember(dta_dims, dim), 1);

% Set the missing values in both datasets to NaN
obs(isnan(mdl)) = NaN;
mdl(isnan(obs)) = NaN;

valid_dta = sum(~isnan(obs), dimpos);

% Compute the number of "valid" data points
% valid_dta = zeros(size(obs));
% valid_dta(~isnan(obs)) = 1;
% valid_dta = sum(valid_dta, dimpos);

for i = 1:length(err_mtrcs)
    switch err_mtrcs{i}
        % 0. Normal errors
        case 'e'
            errs.e  = obs - mdl;
        % 1. Absolute errors
        case 'ae' 
            tmp     = obs - mdl;
            errs.ae = abs(tmp); 
        % 2. Mean absolute errors                                   
        case 'mae'
            tmp      = obs - mdl;
            tmp      = abs(tmp);
            errs.mae = nanmean(tmp, dimpos);

        % 3. Squared errors
        case 'se'
            tmp     = obs - mdl;
            errs.se = tmp.^2;

        % 4. Mean Squared errors
        case 'mse' 
            tmp      = obs - mdl;
            tmp      = tmp.^2;
            errs.mse = nanmean(tmp, dimpos);

        % 5. Root mean squared errors
        case 'rmse'
            tmp       = obs - mdl;
            tmp       = tmp.^2;
            tmp       = nanmean(tmp, dimpos);
            errs.rmse = sqrt(tmp);
        
        % 6. Relative errors
        case 're'
            tmp     = obs - mdl;
            errs.re = tmp./obs;
     
        % 7. Absolute relative errors
        case 'are'
            tmp      = obs - mdl;
            tmp      = tmp./obs;    
            errs.are = abs(tmp);
        
        % 8. Mean absolute relative error (w.r.t. obs)
        case 'mare'
            if find(obs == 0)
                warning('Reference dataset contains zeros!')
            end
            tmp       = abs((obs - mdl)./obs); 
            errs.mare = nanmean(tmp, dimpos);
        
        % 9. Mean absolute relative error (w.r.t. the amplitude of obs)
        case 'mare_alt'
            tmp           = abs(obs - mdl);
            tmp           = nanmean(tmp, dimpos);
            errs.mare_alt = tmp./(max(obs, dimpos) - min(obs, dimpos));    

        % 9. Squared relative error
        case 'sre'
            tmp      = (obs - mdl)./obs;
            errs.sre = tmp.^2;
        
        % 10. Mean squared relative error
        case 'msre'
            if find(obs == 0)
                warning('Reference dataset contains zeros!')
            end
            tmp       = ((obs - mdl)./obs).^2;
            errs.msre = nanmean(tmp, dimpos);
    
        % 11. Root mean squared relative error
        case 'rmsre'
            if find(obs == 0)
                warning('Reference dataset contains zeros!')
            end
            tmp        = ((obs - mdl)./obs).^2
            tmp        = nanmean(tmp, dimpos);
            errs.rmsre = sqrt(tmp);
        
        % 12. Normalized (by the amplitude of obs) root mean square error
        case 'nrmse'
            tmp        = (obs - mdl).^2;
            tmp        = nanmean(tmp, dimpos);
            tmp        = sqrt(tmp);
            errs.nrmse = tmp./(max(obs, dimpos) - min(obs, dimpos)); 
        
        % 13. Coefficient of variation of the RMSE
        case 'cvrmse'
            tmp         = (obs - mdl).^2;
            tmp         = nanmean(tmp, dimpos);
            tmp         = sqrt(tmp);
            errs.cvrmse = tmp./nanmean(obs, dimpos);
        
        % 14. Coefficient of variation of the MAE  
        case 'cvmae'
            tmp        = abs(obs - mdl);
            tmp        = nanmean(tmp, dimpos);
            errs.cvmae = tmp./nanmean(obs, dimpos);    
            
        % 15. Coefficient of determination   
        case 'cod'
            % Compute the mean from both datasets
            mn_obs = nanmean(obs, dimpos);
            mn_mdl = nanmean(mdl, dimpos);
            
            % Remove the mean from the data
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            mdl_cnt = bsxfun(@minus, mdl, mn_mdl);
            
            num   = nansum(obs_cnt.*mdl_cnt, dimpos);
            den   = sqrt(nansum(obs_cnt.^2, dimpos)).*...
                                         sqrt(nansum(mdl_cnt.^2, dimpos));
                                       
            errs.cod = (num./den).^2;
            
        % 16. Percentage bias
        case 'pbias'
            tmp        = obs - mdl;
            errs.pbias = nansum(tmp*100, dimpos)./nansum(obs, dimpos);
        
        % 17. Nash Sutcliffe Efficiency
        case 'nse'
            % Compute the mean of the "observations" (i.e. obs)
            mn_obs  = nanmean(obs, dimpos);
            % Compute the anomalies w.r.t. the mean
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            
            tmp      = obs - mdl;
            a        = nansum(tmp.^2, dimpos);
            b        = nansum(obs_cnt.^2, dimpos);
            
            errs.nse = 1 - a./b;
        
        % 18. Nash Sutcliffe Efficiency (using parameters; e.g. annual cycle)
        case 'nse_param'
            tmp            = obs - mdl;
            a              = nansum(tmp.^2, dimpos);
            b              = nansum((obs - params).^2, dimpos);
            errs.nse_param = 1 - a./b;
        
        % 19. Nash Sutcliffe Efficiency (modified formula (Krause, 2005))
        case 'nse_mod'
            
            mn_obs  = nanmean(obs, dimpos);
            mn_mdl  = nanmean(mdl, dimpos);
        
            obs_cnt = bsxfun(@minus, obs, mn_obs);

            num = nansum(abs(obs - mdl).^params, dimpos);
            den = nansum(obs_cnt.^params, dimpos);
         
            errs.nse_mod = 1 - num./den;
            
        % 19. Nash Sutcliffe Efficiency (alternative formula)
        case 'nse_alt'
            sig_obs = nanstd(obs, dimpos);
            sig_mdl = nanstd(mdl, dimpos);
            
            mn_obs  = nanmean(obs, dimpos);
            mn_mdl  = nanmean(mdl, dimpos);
        
            r = matrixcorr(obs, mdl, dimpos);
            a = sig_mdl./sig_obs;
            b = (mn_mdl - mn_obs)./sig_obs;
         
            errs.nse_alt = 2.*a.*r - a.^2 - b.^2;
            
        % 20. Index of agreement   
        case 'ioa'
            tmp    = obs - mdl;
            % Compute the mean of the observations
            mn_obs = nanmean(obs, dimpos);
            % Compute the anomalies w.r.t. the mean of the OBSERVATIONS
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            mdl_cnt = bsxfun(@minus, mdl, mn_obs);
            
            a      = nansum(tmp.^2, dimpos);
            b      = nansum((abs(mdl_cnt) + abs(obs_cnt)).^2, dimpos);
            
            errs.ioa = 1 - a./b;
            
        % 20. Index of agreement (modified formula (Krause, 2006))
        case 'ioa_mod'
            tmp    = obs - mdl;
            % Compute the mean of the observations
            mn_obs = nanmean(obs, dimpos);
            % Compute the anomalies w.r.t. the mean of the OBSERVATIONS
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            mdl_cnt = bsxfun(@minus, mdl, mn_obs);
            
            a      = nansum(abs(tmp).^params, dimpos);
            b      = nansum((abs(mdl_cnt) + abs(obs_cnt)).^params, dimpos);
            
            errs.ioa_mod = 1 - a./b;
            
        case 'ioa2011'
            mn_obs = nanmean(obs, dimpos);
         
            a = nansum(abs(errs));
            b = 2*nansum(abs(obs - mn_obs));
        
            out = zeros(size(a));
            out(a <= b) = 1 - (a(a <= b))./(b(a <= b));
            out(a >  b) = (b(a > b))./(a(a > b));
        
            errs = out;
            
%     
%     case 'cop'
%         errs = obs - mdl;
%         a    = nansum(errs(2:end, :).^2);
%         b    = nansum((obs(2:end,:) - obs(1:end-1,:)).^2);
%         errs = 1 - a./b;
%         

%         

%         
%     case 'kge_1'
%         
%         r     = nancorr(obs, mdl);
%         mu_s  = nanmean(mdl);
%         mu_o  = nanmean(obs);
%         sig_s = nanstd(mdl);
%         sig_o = nanstd(obs);
%         
%         b     = mu_s./mu_o;
%         g     = (sig_s./mu_s)./(sig_o./mu_o);
%         
%         errs  = 1 - sqrt((r-1).^2 + (b-1).^2 + (g-1).^2);
%         
%         if nargout > 1
%             varargout{1} = r;
%             varargout{2} = b;
%             varargout{3} = g;
%         end
%         
%         
%      case 'kge_2'
%         
%         r     = nancorr(obs, mdl);
%         mu_s  = nanmean(mdl);
%         mu_o  = nanmean(obs);
%         sig_s = nanstd(mdl);
%         sig_o = nanstd(obs);
%         
%         b     = mu_s./mu_o;
%         g     = sig_s./sig_o;
%         
%         errs  = 1 - sqrt((r-1).^2 + (b-1).^2 + (g-1).^2);
%         
%         if nargout > 1
%             varargout{1} = r;
%             varargout{2} = b;
%             varargout{3} = g;
%         end
% 
%     case 'kge_s'
% 
%         r     = nancorr(obs, mdl);
%         mu_s  = nanmean(mdl);
%         mu_o  = nanmean(obs);
%         sig_s = nanstd(mdl);
%         sig_o = nanstd(obs);
%         
%         b     = mu_s./mu_o;
%         g     = sig_s./sig_o;
%         
%         s_r   = params(1);
%         s_b   = params(2);
%         s_g   = params(3);
% 
%         errs = 1 - sqrt((s_r*(r - 1)).^2 + (s_b*(b - 1)).^2 + (s_g*(g - 1)).^2);
%         
%         if nargout > 1
%             varargout{1} = r;
%             varargout{2} = b;
%             varargout{3} = g;
%         end
%            
        % 14. Pearsons correlation coefficient
        case 'corr'
            errs.corr = matrixcorr(obs, mdl, dimpos, true);
            
        % 15. Covariance (normalized by N)   
        case 'cov_0'
            % Compute the mean from both datasets
            mn_obs = nanmean(obs, dimpos);
            mn_mdl = nanmean(mdl, dimpos);
            
            % Remove the mean from the data
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            mdl_cnt = bsxfun(@minus, mdl, mn_mdl);
            
            % Normalize by (N)
            errs.cov = 1./valid_data.*nansum((obs_cnt.*mdl_cnt), dimpos);
            
        % 16. Covariance (normalized by N - 1, BLUE) 
        case 'cov_1'
            % Compute the mean from both datasets
            mn_obs = nanmean(obs, dimpos);
            mn_mdl = nanmean(mdl, dimpos);
            
            % Remove the mean from the data
            obs_cnt = bsxfun(@minus, obs, mn_obs);
            mdl_cnt = bsxfun(@minus, mdl, mn_mdl);
            
            % Nornamlize by (N - 1)
            errs.cov = 1./(valid_data - 1).* ...
                                      nansum((obs_cnt.*mdl_cnt), dimpos);
    end
end

% Throw away the singleton dimension
for i = 1:length(err_mtrcs)
    errs.(err_mtrcs{i}) = squeeze(errs.(err_mtrcs{i}));
end

valid_dta    = squeeze(valid_dta);

varargout{1} = obs_struct;
varargout{2} = mdl_struct;

varargout{3} = tme_out;
varargout{4} = valid_dta;

