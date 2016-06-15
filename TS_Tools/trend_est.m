function [y_out, a, b, sig_p] = trend_est(ts_in, tres_in, indx_row, alpha, wghts, est_method, test_method)
% The function fits a linear trend into each time-series in fld.
% Afterwards, a significance test is performed. The computations follow the
% paper from Santer (2000). 
%--------------------------------------------------------------------------
% Input:        ts_in       [m1 x n]  Matrix which contains the input 
%                                    time-series. 

%                                        

%--------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   March 2015
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
% Updates: 
%--------------------------------------------------------------------------
% References: 
%   B. D. Santer et al. (2000): Statistical significance of
%   trends and trend differences in layer-average atmospheric temperature
%   time series, Journal of Geophysical Research, 105: 7337-7356, 
%   doi: 10.1029/1999JD901105
%   Gilbert, Richard O. (1987), "6.5 Sen's Nonparametric Estimator of
%   Slope", Statistical Methods for Environmental Pollution Monitoring,
%   John Wiley and Sons, pp. 217ñ219, ISBN 978-0-471-28878-7
%--------------------------------------------------------------------------

% Check the input parameters and set default values
if nargin < 7, test_method = 'adjse'; end
if nargin < 6, est_method  = 'TS'; end
if nargin < 5, wghts       = []; end
if nargin < 4, alpha       = 0.95; end
if nargin < 3, indx_row    = true; end
if nargin < 2, tres_in     = 'monthly'; end

% Remove the first row from the time-series, if it contains some indexing
% numbers
[ts_in, tme_vec, num_vec, indx_vec] = separate_data(ts_in, tres_in, indx_row);

% Number of time-steps
nts = size(ts_in, 1);
est_vec = (1:nts)';

if strcmp(est_method, 'LS')
    % Design matrix
    A   = [est_vec ones(nts, 1)];
    
    x_hat = NaN(2, size(ts_in, 2));
    y_hat = NaN(size(ts_in));
    
    for i = 1:size(ts_in, 2)
        % Create an individual A-matrix for each column of observations
        A_tmp = A;
        % Use the ith column of fld as observations
        Y     = ts_in(:, i);
        % Remove NaN-values
        A_tmp(isnan(Y), :) = [];
        Y(isnan(Y))        = [];
        % Save the number of valid observations in a new vector
        nts_nonan(1, i)    = length(Y);
        
        if nts_nonan(1, i) > 0
            % Estimate the parameters through least squares
            x_hat(:, i)     = A_tmp\Y;
            % Compute the trend at all locations in the full A matrix
            y_hat(:, i)     = A*x_hat(:, i);
        end

    end
    
elseif strcmp(est_method, 'TS')
    % Design matrix
    A   = [est_vec ones(nts, 1)];
    mask = zeros(size(ts_in));
    mask(~isnan(ts_in)) = 1;
    nts_nonan = sum(mask);
    
    % Accumulate slopes
    for i = 1:size(ts_in, 1)
        den = repmat(ts_in(i, :), nts, 1) - ts_in;
        num = (est_vec(i) - est_vec)*ones(1, size(ts_in, 2));
       
        C(:, :, i) = den./num;
    end
    
    % Re-arrange the dimensions such that the third dimension corresponds
    % to the different columns of fld
    tmp = permute(C, [1 3 2]);
    % Put all elements of each C(:, :, i)-matrix in a large vector
    C = reshape(tmp, [], size(C,2), 1);        
    
    % Now, the slope can be calculated as the median over a column-vector
    x_hat(1, :) = nanmedian(C, 1);
    
    x_hat(2, :) = nanmedian(ts_in - est_vec*x_hat(1, :));
                       
    y_hat       = A*x_hat;   
    
    x_hat(:, nts_nonan == 0) = NaN;
    y_hat(:, nts_nonan == 0) = NaN;
    
end

% Compute the residuals and their variance
err     = ts_in - y_hat;

if strcmp(test_method, 'naive')
    var_err = sqrt(1./(nts_nonan - 2).*nansum(err.^2));
elseif strcmp(test_method, 'adjse')
    % Compute the lag-1-autocorrelation
    r1 = acf_mtrx(err, 1);
    % Compute the "effective" sample size
    n_e = nts.*(1 - r1)./(1 + r1);
    % Compute the adjusted variance of the residuals
    var_err = sqrt(1./(n_e - 2).*nansum(err.^2));    
end

% Compute the standard error of the trend parameter
std_err   = var_err./sqrt(sum((est_vec - ones(nts, 1)*mean(est_vec)).^2));

% Compute the ratio between the estimated trend and its standard error
t_b = x_hat(1, :)./std_err;

% t_b is distributed as Student's t. Thus, it is compared against a
% critical t-value for a significance level alpha and nts - 2 degrees of
% freedom:
T   = tinv(alpha, nts - 2);

% Significance test: If t_b > T, the null hypothesis b = 0 (i.e. trend is
% not significant) is rejected!
sig_p          = zeros(size(t_b));
sig_p(abs(t_b) > T) = 1;

% Set the elements of sig_p, where no data is available, to NaN
sig_p(:, nts_nonan == 0) = NaN;

a = x_hat(2, :);
b = x_hat(1, :);
    
if strcmp(tres_in, 'annual')
    y_out = [tme_vec y_hat];
else
    y_out = [tme_vec num_vec y_hat];
end

if indx_row == 1
   	indx_vec = [zeros(1, size(y_out, 2) - length(indx_vec)) indx_vec];
    y_out = [indx_vec; y_out];
end










