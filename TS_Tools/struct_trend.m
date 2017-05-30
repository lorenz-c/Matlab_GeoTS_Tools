function [a_out, varargout] = struct_trend(struct_in, vars, alpha, wghts, est_method, test_method, reminsig)
% The function fits a linear trend into each time-series in fld.
% Afterwards, a significance test is performed. The computations follow the
% paper from Santer (2000). 
%--------------------------------------------------------------------------
% Input:        ts_in       [m1 x n]  Matrix which contains the input 
%                                   time-series. 

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
%   John Wiley and Sons, pp. 217?219, ISBN 978-0-471-28878-7
%--------------------------------------------------------------------------

% Check the input parameters and set default values
if nargin < 7, reminsig    = true; end
if nargin < 6, test_method = 'adjse'; end
if nargin < 5, est_method  = 'LS'; end
if nargin < 4, wghts       = []; end
if nargin < 3, alpha       = 0.95; end

% Number of time-steps
nts     = size(struct_in.TimeStamp, 1);
est_vec = (1:nts)';

% Get the data 
fld = struct_in.Data.(vars);

% Get the dimensions of the data
dta_dims = struct_in.Variables.(vars).dimensions;

% Get the "position" of the time-dimension
dimpos = find(ismember(dta_dims, 'time'), 1);

% Compute the numnber of "observations"
nr_dta         = size(fld);
nr_dta(dimpos) = [];

if isempty(nr_dta)
    nr_dta = 1;
end

% Create the "design" matrix 
A  = [est_vec ones(nts, 1)];

% Check if the data is arranged in a 3D-array (e.g. global maps)
if length(dta_dims) == 3
    % If this is the case, we will re-arrange this array into a 2D-matrix.
    bigfld = squeeze(reshape(fld, [], nr_dta(1)*nr_dta(2), 1));
    % Now, the rows correspond to the "time"-dimension (i.e. the
    % observations) and the columns to the variables.
        
    % Calculate the number of real values
    mask              = zeros(size(bigfld));
    mask(~isnan(fld)) = 1;
    nts_nonan         = sum(mask, dimpos);
        
    % If nts_nonan == 0, there might be a pixel with no real values
    % (e.g. over the oceans). We have to remove these columns in
    % advance and save the "position" of the real values for later
    real_indx = find(nts_nonan > 0);
    full_dta  = size(bigfld, 2);
    % Delete the columns with no data
    bigfld(:, nts_nonan == 0) = [];
    nts_nonan(nts_nonan == 0) = [];
    
    fld = bigfld;
end

if strcmp(est_method, 'LS')
    
    % Classical least squares approach
    if dimpos == 2
        fld = fld';
    end
        
    % Calculate number of valid values
    mask              = zeros(size(fld));
    mask(~isnan(fld)) = 1;
    nts_nonan         = sum(mask, 1);
        
%     if max(nts_nonan) == min(nts_nonan)
        % If each variable has the same number of observations, we can
        % calculate the slope and intercept through a sequential
        % estimation.
        % Remove all rows with NaNs from the data, the design matrix,
        % and the time-vector
        est_vec_act = est_vec;
        est_vec_act(isnan(fld(:, 1)), :) = [];
            
        A_act = A;
        A_act(isnan(fld(:, 1)), :) = [];
            
        fld_act = fld;
        fld_act(isnan(fld(:, 1)), :) = [];
            
        x_hat = A_act\fld_act;
%     else
% 
%         % If the variables have not the same number of observations, we
%         % have to loop over each variable and compute the trend
%         % seperately.
%         for i = 1:nr_dta(1)*nr_dta(2)
%             % Create an individual A-matrix for each column of observations
%             A_act = A;
%             % Use the ith column of fld as observations
%             Y     = fld(:, i);      
%             % Remove NaN-values
%             A_act(isnan(Y), :) = [];
%             Y(isnan(Y))        = [];
%                
%             if nts_nonan(1, i) > 0
%                 % Estimate the parameters through least squares
%                 x_hat(:, i)     = A_act\Y;
%             else
%                 x_hat(1:2, i)   = NaN;
%             end          
%         end
%     end
    
    a = x_hat(1, :);
    b = x_hat(2, :);

    y_hat             = A*x_hat;
    y_hat(isnan(fld)) = NaN;  
    
elseif strcmp(est_method, 'TS')
  
    % Calculate number of valid values
    mask              = zeros(size(fld));
    mask(~isnan(fld)) = 1;
    nts_nonan         = sum(mask, dimpos);
    
    % If the columns correspond to the time-dimension, transpose the
    % data matrix
    if dimpos == 2
        fld = fld';
    end
       
    % Use two different approaches, depending on the number of
    % datapoints and the number of time steps. There might be no issues
    % for small amounts of data. But if e.g. the number of variables
    % gets quite large (say, > 10000), we might run into memory issues.
    % Therefore, the two approaches might seem a bit counter-intuitive
    % as we're looping over the larger dimensions. But that saves
    % actually a lot of memory. 
    if nr_dta >= max(nts_nonan)
        % If there are more variables than observations, then we'll
        % loop over the number of variables and throw away all missing
        % data. Then, we'll calculate the slope for each variable
        % within the loop
        for i = 1:size(fld, 2)
            % First, truncate the data and the time vector so that it 
            % only contais valid values
            fld_act     = fld(:, i);
            est_vec_act = est_vec;
            
            est_vec_act(isnan(fld_act)) = [];
            fld_act(isnan(fld_act))     = [];
            nts_act                     = length(fld_act);
                
            % We then generate two long vectors for both the
            % observations and the time values.
            fld_1 = repmat(fld_act, nts_act, 1);
            fld_2 = repmat(fld_act', nts_act, 1);
            num   = fld_1 - fld_2(:);
            
            est_1 = repmat(est_vec_act, nts_act, 1);
            est_2 = repmat(est_vec_act', nts_act, 1);
            den   = est_1 - est_2(:);
            
            % The full ensemble of slopes is than generated by dividing
            % the two differences 
            C = num./den;
            % Finally, compute the slope from the median of the
            % ensemble
            a(i) = nanmedian(C);
        end
        
    else
        % If there are more observations than variables, we'll loop
        % over the number of observations and compute (for every
        % time-step) a 3D-array. This array, however, might get very
        % large if there are many observations and many variables.
        for i = 1:nts
            C(:, :, i) = bsxfun(@rdivide, bsxfun(@minus, fld(i, :), ...
                    fld), est_vec(i) - est_vec);
        end
        
        % Then, we'll permute the C-matrix, so that the first and
        % second dimension correspond to the ensemble of slopes
        Cprm = permute(C, [1 3 2]);
        % Re-arrage the 3D-array into a 2D-matrix
        Cprm = reshape(Cprm, [],size(C,2),1 );
        % Compute the median of slopes from that matrix
        a = nanmedian(Cprm);        
    end
       
    % Compute the intercept (Helsel and Hirsch 1995, page 267; 
    % Conover 1999, page 336)
    b = nanmedian(fld) - a.*nanmedian(est_vec);
        
    % Compute the y-values of the trend line
    y_hat = A*[a; b];   
    
    % Set the parameters where no observations are available to NaN
    a(nts_nonan == 0) = NaN;
    b(nts_nonan == 0) = NaN;
        
    y_hat(isnan(fld)) = NaN;
end

% Compute the residuals and their variance
err = fld - y_hat;

if strcmp(test_method, 'naive')
    var_err = sqrt(1./(nts_nonan - 2).*nansum(err.^2));
    
elseif strcmp(test_method, 'adjse')
    % Compute the lag-1-autocorrelation
    r1 = acf_mtrx(err, 1);
    
    % Compute the "effective" sample size
    n_e = nts_nonan.*(1 - r1)./(1 + r1);
    
    % Compute the adjusted variance of the residuals
    var_err = sqrt(1./(n_e - 2).*nansum(err.^2));    
end

% Compute the standard error of the trend parameter
std_err   = var_err./sqrt(sum((est_vec - ones(nts, 1)*mean(est_vec)).^2));

% Compute the ratio between the estimated trend and its standard error
t_b = a./std_err;

% t_b is distributed as Student's t. Thus, it is compared against a
% critical t-value for a significance level alpha and nts - 2 degrees of
% freedom:
T   = tinv(alpha, nts - 2);

% Significance test: If t_b > T, the null hypothesis b = 0 (i.e. trend is
% not significant) is rejected!
sig_p               = zeros(size(t_b));
sig_p(abs(t_b) > T) = 1;

% Set the elements of sig_p, where no data is available, to NaN
sig_p(:, nts_nonan == 0) = NaN;


if length(dta_dims) == 3
    % If the original data had 3 dimensions (--> maps), we have to
    % reconstruct these maps 
    tmp               = NaN(nts, nr_dta(1)*nr_dta(2));
    tmp(:, real_indx) = y_hat;
    y_hat             = reshape(tmp, nts, nr_dta(1), nr_dta(2));
    
    a_out            = NaN(nr_dta(1), nr_dta(2));
    a_out(real_indx) = a;
    
    if strcmp(est_method, 'LS')
        b_out = NaN(nr_dta(1), nr_dta(2));
        b_out(real_indx) = b;
    end
    
    sig_out            = NaN(nr_dta(1), nr_dta(2));
    sig_out(real_indx) = sig_p;
else
    a_out   = a;
    sig_out = sig_p;
    
    if strcmp(est_method, 'LS'), b_out = b; end
end



if reminsig == true
    a_out     = a_out  .* sig_out;
    
    if strcmp(est_method, 'LS')
        b_out     = b_out  .* sig_out;
    end
end

ts_out = y_hat;
    



varargout{1} = sig_out;
varargout{2} = ts_out;
if strcmp(est_method, 'LS')
    varargout{3} = b_out;
end
    




