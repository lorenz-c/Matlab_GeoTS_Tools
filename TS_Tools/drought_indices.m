function out = drought_indices(inpt, scale, type)

if strcmp(type, 'SPI_daily')
    
elseif strcmp(type, 'SPI_param_monthly')
    % Fit a gamma distribution to the input data
    pd = fitdist(inpt, 'Gamma');
elseif strcmp(type, 'SPI_emp_monthly')
    tmp = NaN(length(inpt) - scale + 1, scale);
    
    for i = 1:scale
        tmp(:, i) = inpt(i:end-scale + i);
    end
    
    % Compute seasonal sums
    ssnl_sm = sum(tmp, 2);
    
    
    
end