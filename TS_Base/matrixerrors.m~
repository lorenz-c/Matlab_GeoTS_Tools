function [errs, varargout] = matrixerrors(fld1, quant, varargin)
%--------------------------------------------------------------------------
% The function is similar to the struct_errors function, but does not
% require any metadata or other ancillary information. It therefore does
% not perform any subsetting of the data. 
% The first input (fld1) is assumed to consist e.g. some observations,
% against which other input datasets (in varargin) are compared. 
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

if nargin < 3, quant    = 'rmse';    end


for i = 1:length(varargin)
    if size(fld1) ~= size(varargin{i})
        error('Datasets must have the same dimensions!')  
    else
        fld_val{i} = varargin{i};
        fld_val{i}(isnan(fld1)) = NaN;
        
        fld1(isnan(fld_val{i})) = NaN;
    end
end
    
nts  = length(fld1(:,1));

for i = 1:length(varargin)
    switch quant
    
        % 1. Absolute errors
        case 'ae' 
            errs = fld1 - fld2;
            errs = abs(errs); 

        % 2. Mean absolute errors                                   
        case 'mae'
            errs = fld1 - fld2;
            errs = abs(errs);
            errs = nanmean(errs);

        % 3. Squared errors
        case 'se'
            errs = fld1 - fld2;
            errs = errs.^2;

        % 4. Mean Squared errors
        case 'mse' 
            errs = fld1 - fld2;
            errs = errs.^2;
            errs = nanmean(errs);

        % 5. Root mean squared errors
        case 'rmse'
            errs = fld1 - fld2;
            errs = errs.^2;
            errs = nanmean(errs);
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
        
        % 8. Mean absolute relative errors
        case 'mare'
            errs = fld1 - fld2;
            fld1(fld1 == 0) = NaN;
            errs = abs(errs);
            errs = errs./fld1;    
            errs = nanmean(errs);
        
        case 'mare_2'
            errs = fld1 - fld2;
            errs = abs(errs);
            errs = nanmean(errs);
            errs = errs./(max(fld1) - min(fld1));    

        % 9. Squared relative errors
        case 'sre'
            errs = fld1 - fld2;
            errs = errs./fld1;
            errs = errs.^2;
        
        % 10. Mean squared relative errors
        case 'msre'
            errs = fld1 - fld2;
            errs = errs./fld1;
            errs = errs.^2;
            errs = nanmean(errs);
    
        % 11. Root mean squared relative erros
        case 'rmsre'
            errs = fld1 - fld2;
            errs = errs./fld1;
            errs = errs.^2;
            errs = nanmean(errs);
            errs = sqrt(errs);
        
        % 12. Normalized root mean square deviation
        case 'nrmse'
            errs = fld1 - fld2;
            errs = errs.^2;
            errs = nanmean(errs);
            errs = sqrt(errs);
            errs = errs./(max(fld1) - min(fld1)); 
        
        % 13. Coefficient of variation of the RMSE
    	case 'cvrmse'
            errs = fld1 - fld2;
            errs = errs.^2;
            errs = nanmean(errs);
            errs = sqrt(errs);
            errs = errs./nanmean(fld1);
        
        case 'cvmae'
            errs = fld1 - fld2;
            errs = abs(errs);
            errs = nanmean(errs);
            errs = errs./nanmean(fld1);    
        
        case 'corr'
            errs = nancorr(fld1, fld2);
        
        case 'rsq'
            errs = (nancorr(fld1, fld2)).^2;
        
        case 'cod'
            o_anom = fld1 - ones(nts, 1)*nanmean(fld1);
            p_anom = fld2 - ones(nts, 1)*nanmean(fld2);
        
            nom   = nansum(o_anom.*p_anom);
            denom = sqrt( nansum(o_anom.^2)).*sqrt(nansum(p_anom.^2));
        
            errs  = (nom./denom).^2;
        
        case 'pbias'
            errs = fld1 - fld2;
            errs = nansum(errs*100)./nansum(fld1);
    
        case 'pbias2'
            errs = fld2 - fld1;
            errs = nansum(errs*100)./nansum(fld1);
        
        case 'nse'
            errs = fld1 - fld2;
            a    = nansum(errs.^2);
            b    = nansum((fld1 - ones(nts,1)*nanmean(fld1)).^2);
            errs = 1 - a./b;
        
        case 'nse_m'
            errs = fld1 - fld2;
            a    = nansum(errs.^2);
        
            [tmp, mnth_mn]       = remsc(fld1, [], 0, 3);
            nr_yrs               = size(fld1, 1)/12;
            mnth_mn              = repmat(mnth_mn, nr_yrs, 1);
            mnth_mn(isnan(fld1)) = NaN;

            b = nansum((fld1 - mnth_mn).^2);
        
            errs = 1 - a./b;
      
        case 'nse_alt'
            errs = fld1 - fld2;
            a    = nansum(errs.^2);
            b    = nansum((fld1 - params).^2);
        
            errs = 1 - a./b;
        
        case 'nse2'
        
            sig_s = nanstd(fld2);
            sig_o = nanstd(fld1);
        
            mu_s  = nanmean(fld2);
            mu_o  = nanmean(fld1);
        
            r = nancorr(fld1, fld2);
            a = sig_s./sig_o;
            b = (mu_s - mu_o)./sig_o;
        
            errs = 2.*a.*r - a.^2 - b.^2;
    
        case 'cop'
            errs = fld1 - fld2;
            a    = nansum(errs(2:end, :).^2);
            b    = nansum((fld1(2:end,:) - fld1(1:end-1,:)).^2);
            errs = 1 - a./b;
        
        case 'ioa'
            errs   = fld1 - fld2;
            mn_obs = ones(nts,1)*nanmean(fld1);
            a      = nansum(errs.^2);
            b      = nansum((abs(fld2 - mn_obs) + abs(fld1 - mn_obs)).^2);
            errs   = 1 - a./b;
        
        case 'ioa2011'
            mn_obs = ones(nts,1)*nanmean(fld1);
        
            a = nansum(abs(errs));
            b = 2*nansum(abs(fld1 - mn_obs));
        
            out = zeros(size(a));
            out(a <= b) = 1 - (a(a <= b))./(b(a <= b));
            out(a >  b) = (b(a > b))./(a(a > b));
        
            errs = out;
        
        case 'kge_1'
        
            r     = nancorr(fld1, fld2);
            mu_s  = nanmean(fld2);
            mu_o  = nanmean(fld1);
            sig_s = nanstd(fld2);
            sig_o = nanstd(fld1);
        
            b     = mu_s./mu_o;
            g     = (sig_s./mu_s)./(sig_o./mu_o);
        
            errs  = 1 - sqrt((r-1).^2 + (b-1).^2 + (g-1).^2);
        
            if nargout > 1
                varargout{1} = r;
                varargout{2} = b;
                varargout{3} = g;
            end
        
        
        case 'kge_2'
        
            r     = nancorr(fld1, fld2);
            mu_s  = nanmean(fld2);
            mu_o  = nanmean(fld1);
            sig_s = nanstd(fld2);
            sig_o = nanstd(fld1);
        
            b     = mu_s./mu_o;
            g     = sig_s./sig_o;
        
            errs  = 1 - sqrt((r-1).^2 + (b-1).^2 + (g-1).^2);
        
            if nargout > 1
                varargout{1} = r;
                varargout{2} = b;
                varargout{3} = g;
            end

        case 'kge_s'

            r     = nancorr(fld1, fld2);
            mu_s  = nanmean(fld2);
            mu_o  = nanmean(fld1);
            sig_s = nanstd(fld2);
            sig_o = nanstd(fld1);
        
            b     = mu_s./mu_o;
            g     = sig_s./sig_o;
        
            s_r   = params(1);
            s_b   = params(2);
            s_g   = params(3);

            errs = 1 - sqrt((s_r*(r - 1)).^2 + (s_b*(b - 1)).^2 + ...
                                                         (s_g*(g - 1)).^2);
        
            if nargout > 1
                varargout{1} = r;
                varargout{2} = b;
                varargout{3} = g;
            end
    end
    
end




        
