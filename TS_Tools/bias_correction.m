function obs_p = bias_correction(mdl_h, obs_h, mdl_p, method, cdf_method, thresh, varargin)
  


% Functions implemented from Cannon et al. (2015): Bias Correction of GCM 
% Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes 
% in Quantiles and Extremes?, doi: 10.1175/JCLI-D-14-00754.1

% For e.g. daily precipitation with it's mixed discrete-continous nature,
% small values should be removed or replaced. Therefore, the function
% replaces the values which are smaller than thresh with uniform random
% values between zero and the threshold.
if thresh > 0
    zeros_mdl_h = find(mdl_h < thresh);
    zeros_mdl_p = find(mdl_p < thresh);
    zeros_obs_h = find(obs_h < thresh);


    rand_mdl_h  = thresh*rand(length(zeros_mdl_h), 1);
    rand_mdl_p  = thresh*rand(length(zeros_mdl_p), 1);
    rand_obs_h  = thresh*rand(length(zeros_obs_h), 1);

    mdl_h(zeros_mdl_h) = rand_mdl_h;
    mdl_p(zeros_mdl_p) = rand_mdl_p;
    obs_h(zeros_obs_h) = rand_obs_h;
end




if strcmp(method, 'qm') 
    
   
    [y_mdl_h, x_mdl_h] = cdf_transform(mdl_h, cdf_method, varargin);
    [y_obs_h, x_obs_h] = cdf_transform(obs_h, cdf_method, varargin);
     
    % Get the probabilities of the values during the validation period
    y_mdl_p = interp1(x_mdl_h, y_mdl_h, mdl_p);
    
    % 
    obs_p(:, 1) = quantile(obs_h, y_mdl_p);
   
elseif strcmp(method, 'dqm')
    
    a = (nanmean(mdl_h)*mdl_p)/nanmean(mdl_p);
    b = nanmean(mdl_p)/nanmean(mdl_h);
    
    [y_mdl_h, x_mdl_h] = cdf_transform(mdl_h, cdf_method, varargin);
   
    y_mdl_p = interp1(x_mdl_h, y_mdl_h, a);

    obs_p = quantile(obs_h, y_mdl_p)*b;
    
elseif strcmp(method, 'qdm')
    
    [y_mdl_p, x_mdl_p] = cdf_transform(mdl_p, cdf_method, varargin);
    y_tau_p            = interp1(x_mdl_p, y_mdl_p, mdl_p);
    delta_m            = mdl_p./quantile(mdl_h, y_tau_p);  
    x_omhp             = quantile(obs_h, y_tau_p);    
    obs_p              = x_omhp.*delta_m;
  
end

if thresh > 0
    obs_p(zeros_mdl_p) = 0;
end

end



function [FX, X] = cdf_transform(x, method, varargin)

    if strcmp(method, 'ecdf')
        [FX, X] = ecdf(x);
        FX      = FX(2:end);
        X       = X(2:end);
    elseif strcmp(cdf_method, 'ksdensity')
        [FX, X] = ksdensity(x, 'function', 'cdf', ...
                                                    'Support', 'positive');
    elseif strcmp(cdf_method, 'pareto')
        dst = paretotails(x, varargin{1}, varargin{2});
        X   = linspace(min(x), max(x), 100);
        FX  = dst.cdf(x);  
    end
 end
    
 
