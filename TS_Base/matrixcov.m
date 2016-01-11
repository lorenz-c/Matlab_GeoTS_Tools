function C = matrixcov(X, Y, method, div, remmn, nomiss, corrneg)
% nancov_matrix computes the column-wise covariance matrix between X and Y.
% NaNs are treated as missing elements. The user can chose between two
% methods:
% - pairwise: for each pair of columns, the NaNs are removed separately,
%             which means that the covariances are based on different 
%             (e.g. time) periods
% - row:      if a row in X or Y contains NaN-values, it is removed
%             completely. Thus, the covariances are based on a consistent 
%             period, but some data might be unused
% Furthermore, the covariances are either normalized by N (div = 0) or 
% N - 1 (div = 1, default) to compute either the second moment matrix of 
% the observations about ther mean (div = 0) or the best unbiased estimate 
% of the covariance matrix (div = 1) if the observations are from a normal
% distribution.

if nargin < 7, corrneg = false; end
if nargin < 6, nomiss = 0; end
if nargin < 5, remmn = 1; end
if nargin < 4, div = 0; end
if nargin < 3, method = 'pairwise'; end

% Check the dimensions of the input matrices
[r1, c1] = size(X);
[r2, c2] = size(Y);

if r1 ~= r2 | c1 ~= c2
    error('The input matrices must be of the same size')
end

if nomiss == 0
    if strcmp(method, 'pairwise')
        % Set all missing elements in either X or Y to NaN
        X(isnan(Y)) = NaN;
        Y(isnan(X)) = NaN;
        
        if remmn == 1
            % Compute the mean from both X and Y
            mn_X = nanmean(X);
            mn_Y = nanmean(Y);
            
            % Remove the mean from both X and Y
            X    = X - repmat(mn_X, r1, 1);
            Y    = Y - repmat(mn_Y, r1, 1);
        end
       
        % Compute a mask, which is used to count the number of valid values
        mask              = zeros(size(X));
        mask(~isnan(X)) = 1;
         
        div_mat         = mask'*mask;
        div_mat         = div_mat - div;
       
        % Set all NaN-elements to 0
        X(isnan(X)) = 0;
        Y(isnan(Y)) = 0;

        C = (1./div_mat).*(X'*Y);

    elseif strcmp(method, 'row')
        % Mask contains ones at the positions where either X or Y has
        % NaN-elements.
        mask           = zeros(size(X));
        mask(isnan(X)) = 1;
        mask(isnan(Y)) = 1;
    
        % Compute the row-wise sum of mask. The rows where sum(mask, 2) > 0 are
        % removed from both X and Y.
        mask = sum(mask, 2);
        X(mask > 0, :) = [];
        Y(mask > 0, :) = [];
    
        % Compute the number of observations after removing the NaN-elements
        N = size(X, 1);
    
        if remmn == 1
            % Compute the mean of X and Y
            mn_X = mean(X, 1);
            mn_Y = mean(Y, 1);
    
            % Remove the mean from the input matrices
            X = X - repmat(mn_X, N, 1);
            Y = Y - repmat(mn_Y, N, 1);
        end
    
        % Compute the empirical covariance matrix between X and Y
        if div == 1
            C = (1/(N-1))*(X'*Y);
        elseif div == 0
            C = (1/(N))*(X'*Y);
        end
    end
    
elseif nomiss == 1
    
    % Compute the number of observations after removing the NaN-elements
    N = size(X, 1);
    
    if remmn == 1
        % Compute the mean of X and Y
        mn_X = mean(X, 1);
        mn_Y = mean(Y, 1);
     
        % Remove the mean from the input matrices
        X = X - repmat(mn_X, N, 1);
        Y = Y - repmat(mn_Y, N, 1);
    end

    % Compute the empirical covariance matrix between X and Y
	if div == 1
        C = (X'*Y)/(N-1);
	elseif div == 0
        C = (X'*Y)/N;
    end
end



if corrneg == true
    
    EPS = 10^-6;
    ZERO = 10^-10;
    
    [V, D]       = eig(C);
    
    if min(D) < ZERO
        D(D <= ZERO) = EPS;
        C            = V*diag(diag(C))*V';
    end
    
end

    
    
