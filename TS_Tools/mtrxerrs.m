function [errs, varargout] = mtrxerrs(fld1, fld2, quant, tscale, params, cntr)


if nargin < 6, cntr     = 0;         end
if nargin < 5, params   = [1 1 1];   end
if nargin < 4, tscale   = 'monthly'; end
if nargin < 3, quant    = 'rmse';    end

if strcmp(tscale, 'monthly')
    dta_clms = [2 4];
elseif strcmp(tscale, 'daily')
    dta_clms = [2 5];
elseif strcmp(tscale, 'none')
    dta_clms = [1 1];
end

if size(fld1, 1) ~= size(fld2, 1)
    fld2 = find_ref_dtes(fld1, tscale, fld2);
end
    


fld1 = fld1(dta_clms(1):end, dta_clms(2):end);
fld2 = fld2(dta_clms(1):end, dta_clms(2):end);



fld1(isnan(fld2)) = NaN;
fld2(isnan(fld1)) = NaN;


nts  = length(fld1(:,1));

if cntr == 1
    fld1 = fld1 - ones(nts, 1)*nanmean(fld1);
    fld2 = fld2 - ones(nts, 1)*nanmean(fld2);
end

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

        errs = 1 - sqrt((s_r*(r - 1)).^2 + (s_b*(b - 1)).^2 + (s_g*(g - 1)).^2);
        
        if nargout > 1
            varargout{1} = r;
            varargout{2} = b;
            varargout{3} = g;
        end
               
end




        
