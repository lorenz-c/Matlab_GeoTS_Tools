function error_nme = full_error_names(short_names)

if isstr(short_names)
    nms{1} = short_names;
else
    nms    = short_names;
end


for i = 1:length(nms)
    if strcmp(nms{i}, 'ae')
        error_nme{i} = 'Absolute error';
    elseif strcmp(nms{i}, 'mae')
        error_nme{i} = 'Mean absolute error';
    elseif strcmp(nms{i}, 'se')
        error_nme{i} = 'Squared error';
    elseif strcmp(nms{i},'mse')
        error_nme{i} = 'Mean squared error';
    elseif strcmp(nms{i}, 'rmse')
        error_nme{i} = 'Root mean squared error';
    elseif strcmp(nms{i}, 're')
        error_nme{i} = 'Relative errors';
    elseif strcmp(nms{i}, 'are')
        error_nme{i} = 'Absolute relative errors';
    elseif strcmp(nms{i}, 'mare')
        error_nme{i} = 'Mean absolute relative errors';
    elseif strcmp(nms{i}, 'mare_alt')
        error_nme{i} = 'Mean absolute relative errors (w.r.t. amplitude)';
    elseif strcmp(nms{i}, 'sre')
        error_nme{i} = 'Squared relative errors';
    elseif strcmp(nms{i}, 'msre')
        error_nme{i} = 'Mean squared relative errors';
    elseif strcmp(nms{i}, 'rmsre')
        error_nme{i} = 'Root mean squared relative errors';    
    elseif strcmp(nms{i}, 'nrmse')
        error_nme{i} = ...
                  'Normalized root mean squared errors (w.r.t. amplitude)';
    elseif strcmp(nms{i}, 'cvrmse')
        error_nme{i} = 'Coefficient of variation of the RMSE';        
    elseif strcmp(nms{i}, 'cvmae')
        error_nme{i} = 'Coefficient of variation of the MAE';   
    elseif strcmp(nms{i}, 'cod') 
        error_nme{i} = 'Coefficient of determination';
    elseif strcmp(nms{i}, 'pbias')
        error_nme{i} = 'Percentage bias';
    elseif strcmp(nms{i}, 'nse')
        error_nme{i} = 'Nash Sutcliffe Efficiency';
    elseif strcmp(nms{i}, 'nse_param')
        error_nme{i} = 'Nash Sutcliffe Efficiency w.r.t. parameter';
    elseif strcmp(nms{i}, 'nse_mod')
        error_nme{i} = 'Modified NSE (Krause, 2006)';
    elseif strcmp(nms{i}, 'nse_alt')
        error_nme{i} = 'Alternative NSE';
    elseif strcmp(nms{i}, 'ioa')
        error_nme{i} = 'Index of agreement';
    elseif strcmp(nms{i}, 'ioa_mod')
        error_nme{i} = 'Modified index of agreement (Krause, 2006)';
    elseif strcmp(nms{i}, 'ioa2011')
        error_nme{i} = 'Modified index of agreement (Willmott, 2011)';
    elseif strcmp(nms{i}, 'corr')
        error_nme{i} = 'Pearsons correlation coefficient';
    elseif strcmp(nms{i}, 'cov_0')
        error_nme{i} = 'Covariance, normalized by N';
    elseif strcmp(nms{i}, 'cov_1')
        error_nme{i} = 'Covariance, normalized by N - 1';
    end
end