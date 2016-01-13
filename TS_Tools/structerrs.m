function [errs, varargout] = structerrs(struct1, struct2, vars, quant, dim, params)
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
% - pbias     Relative bias (w.r.t. struct1)
%
% Depending on the chosen error metric, the output will have the same
% dimensions as the input data. If 'average' error metrics are required,
% the dimension along which the errors are computed is removed.
%
% The function already contains some placeholders for more metrics, which
% will be successively added.
%--------------------------------------------------------------------------
% INPUT:
% - struct1   CF datastructure (e.g. with observations)
% - struct2   CF datastructure
% - vars      Variable for which the errors should be computed 
%             (e.g. 'prec')
% - quant     Error metric (see above)
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

if nargin < 7, cntr     = 0;         end
if nargin < 6, params   = [1 1 1];   end
if nargin < 5, dim      = 'time';    end
if nargin < 4, quant    = 'rmse';    end


% Truncate the two structures to the same time period
strt_yr = max([struct1.Data.time(1, 1), struct2.Data.time(1, 1)]);
end_yr  = min([struct1.Data.time(end, 1), struct2.Data.time(end, 1)]);

struct1 = trunc_TS(struct1, strt_yr, end_yr);
struct2 = trunc_TS(struct2, strt_yr, end_yr);

% Save the reference period in tme_out
tme_out.Data      = struct1.Data.time;
tme_out.TimeStamp = struct1.TimeStamp;

% Get the data for the desired variable from both datasets
fld1 = struct1.Data.(vars);
fld2 = struct2.Data.(vars);

% Get the dimensions of the data
dta_dims = struct1.Variables.(vars).dimensions;

% Get the "position" of the dimension, over which the errors should be
% calculated
dimpos = find(ismember(dta_dims, dim), 1);

% Get the data
fld1 = struct1.Data.(vars);
fld2 = struct2.Data.(vars);

% Set the missing values in both datasets to NaN
fld1(isnan(fld2)) = NaN;
fld2(isnan(fld1)) = NaN;

% Compute the number of "valid" data points
valid_dta = zeros(size(fld1));
valid_dta(~isnan(fld1)) = 1;
valid_dta = sum(valid_dta, dimpos);


switch quant
    
    % 1. Absolute errors
    case 'ae' 
        errs = fld1 - fld2;
        errs = abs(errs); 

    % 2. Mean absolute errors                                   
    case 'mae'
        errs = fld1 - fld2;
        errs = abs(errs);
        errs = nanmean(errs, dimpos);

    % 3. Squared errors
    case 'se'
        errs = fld1 - fld2;
        errs = errs.^2;

    % 4. Mean Squared errors
    case 'mse' 
        errs = fld1 - fld2;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);

    % 5. Root mean squared errors
    case 'rmse'
        errs = fld1 - fld2;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);
        errs = sqrt(errs);
        
    % 6. Relative errors
    case 're'
        errs = fld1 - fld2;
        errs = errs./fld1;
     
    % 7. Absolute relative errors
    case 'are'
        errs = fld1 - fld2;
        errs = errs./fld1;    
        errs = abs(errs);
        
    % 8. Mean absolute relative error (w.r.t. fld1)
    case 'mare'
        errs = fld1 - fld2;
        fld1(fld1 == 0) = NaN;
        errs = abs(errs);
        errs = errs./fld1;    
        errs = nanmean(errs, dimpos);
        
    % 9. Mean absolute relative error (w.r.t. the amplitude of fld1)
    case 'mare_2'
        errs = fld1 - fld2;
        errs = abs(errs);
        errs = nanmean(errs, dimpos);
        errs = errs./(max(fld1, dimpos) - min(fld1, dimpos));    

    % 9. Squared relative error
    case 'sre'
        errs = fld1 - fld2;
        errs = errs./fld1;
        errs = errs.^2;
        
    % 10. Mean squared relative error
    case 'msre'
        errs = fld1 - fld2;
        errs = errs./fld1;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);
    
    % 11. Root mean squared relative error
    case 'rmsre'
        errs = fld1 - fld2;
        errs = errs./fld1;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);
        errs = sqrt(errs);
        
    % 12. Normalized (by the amplitude of fld1) root mean square error
    case 'nrmse'
        errs = fld1 - fld2;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);
        errs = sqrt(errs);
        errs = errs./(max(fld1, dimpos) - min(fld1, dimpos)); 
        
    % 13. Coefficient of variation of the RMSE
    case 'cvrmse'
        errs = fld1 - fld2;
        errs = errs.^2;
        errs = nanmean(errs, dimpos);
        errs = sqrt(errs);
        errs = errs./nanmean(fld1, dimpos);
        
    % 13. Coefficient of variation of the MAE  
    case 'cvmae'
        errs = fld1 - fld2;
        errs = abs(errs);
        errs = nanmean(errs, dimpos);
        errs = errs./nanmean(fld1, dimpos);    
        
%     case 'corr'
%         errs = nancorr(fld1, fld2);
%         
%     case 'rsq'
%         errs = (nancorr(fld1, fld2)).^2;
%         
%     case 'cod'
%         o_anom = fld1 - ones(nts, 1)*nanmean(fld1);
%         p_anom = fld2 - ones(nts, 1)*nanmean(fld2);
%         
%         nom   = nansum(o_anom.*p_anom);
%         denom = sqrt( nansum(o_anom.^2)).*sqrt(nansum(p_anom.^2));
%         
%         errs  = (nom./denom).^2;
        
    case 'pbias'
        errs = fld1 - fld2;
        errs = nansum(errs*100, dimpos)./nansum(fld1, dimpos);
           
%     case 'nse'
%         errs = fld1 - fld2;
%         a    = nansum(errs.^2, dimpos);
%         b    = nansum((fld1 - ones(nts,1)*nanmean(fld1)).^2);
%         errs = 1 - a./b;
        
%     case 'nse_m'
%         errs = fld1 - fld2;
%         a    = nansum(errs.^2);
%         
%         [tmp, mnth_mn]       = remsc(fld1, [], 0, 3);
%         nr_yrs               = size(fld1, 1)/12;
%         mnth_mn              = repmat(mnth_mn, nr_yrs, 1);
%         mnth_mn(isnan(fld1)) = NaN;
% 
%         b = nansum((fld1 - mnth_mn).^2);
%         
%         errs = 1 - a./b;
%       
%     case 'nse_alt'
%         errs = fld1 - fld2;
%         a    = nansum(errs.^2);
%         b    = nansum((fld1 - params).^2);
%         
%         errs = 1 - a./b;
        
%     case 'nse2'
%         
%         sig_s = nanstd(fld2, dimpos);
%         sig_o = nanstd(fld1, dimpos);
%         
%         mu_s  = nanmean(fld2, dimpos);
%         mu_o  = nanmean(fld1, dimpos);
%         
%         r = nancorr(fld1, fld2);
%         a = sig_s./sig_o;
%         b = (mu_s - mu_o)./sig_o;
%         
%         errs = 2.*a.*r - a.^2 - b.^2;
%     
%     case 'cop'
%         errs = fld1 - fld2;
%         a    = nansum(errs(2:end, :).^2);
%         b    = nansum((fld1(2:end,:) - fld1(1:end-1,:)).^2);
%         errs = 1 - a./b;
%         
%     case 'ioa'
%         errs   = fld1 - fld2;
%         mn_obs = ones(nts,1)*nanmean(fld1);
%         a      = nansum(errs.^2);
%         b      = nansum((abs(fld2 - mn_obs) + abs(fld1 - mn_obs)).^2);
%         errs   = 1 - a./b;
%         
%     case 'ioa2011'
%         mn_obs = ones(nts,1)*nanmean(fld1);
%         
%         a = nansum(abs(errs));
%         b = 2*nansum(abs(fld1 - mn_obs));
%         
%         out = zeros(size(a));
%         out(a <= b) = 1 - (a(a <= b))./(b(a <= b));
%         out(a >  b) = (b(a > b))./(a(a > b));
%         
%         errs = out;
%         
%     case 'kge_1'
%         
%         r     = nancorr(fld1, fld2);
%         mu_s  = nanmean(fld2);
%         mu_o  = nanmean(fld1);
%         sig_s = nanstd(fld2);
%         sig_o = nanstd(fld1);
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
%         r     = nancorr(fld1, fld2);
%         mu_s  = nanmean(fld2);
%         mu_o  = nanmean(fld1);
%         sig_s = nanstd(fld2);
%         sig_o = nanstd(fld1);
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
%         r     = nancorr(fld1, fld2);
%         mu_s  = nanmean(fld2);
%         mu_o  = nanmean(fld1);
%         sig_s = nanstd(fld2);
%         sig_o = nanstd(fld1);
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
end

% Throw away the singleton dimensions
errs      = squeeze(errs);
valid_dta = squeeze(valid_dta);

varargout{1} = tme_out;
varargout{2} = valid_dta;
