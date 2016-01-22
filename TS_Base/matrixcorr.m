function R = matrixcorr(M1, M2, dim, remnan);
%--------------------------------------------------------------------------
% The function computes the column- or row-wise correlation between two
% matrices. 
%--------------------------------------------------------------------------
% INPUT:
% - M1, M2    Matrices which must have the same dimension
% - dim       Dimension from which the correlation should be calculated
%             (can be set to 1 (default) or 2).
% - remnan    Boolean parameter. If set to true (default), the function 
%             removes the NaNs from the matrices. If set to false, the 
%             function assumes that there are no missing values in the 
%             matrices. Otherwise, the function gives NaNs in colums 
%             (or rows) where one or more values are missing.
%--------------------------------------------------------------------------
% OUTPUT:
% - s         Row- or column-vector with correlations
%--------------------------------------------------------------------------
% EXAMPLE:
% >> R = matrixcorr(M1, M2, 2, true)
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

if nargin < 4, remnan = true; end
if nargin < 3, dim = 1; end

% Check if the matrices have the same size
if size(M1) ~= size(M2)
    error('Matrices must have the same size!')
end

if remnan == true
    % Check for missing values in the matrices
    M1(isnan(M2)) = NaN;
    M2(isnan(M1)) = NaN;
    % Compute the nanmean
    mn1  = nanmean(M1, dim);
    mn2  = nanmean(M2, dim);    
    % Set all NaN-elements to zero for the following matrix multiplication
	M1(isnan(M1)) = 0;
    M2(isnan(M2)) = 0;   
else
    % Compute the mean
    mn1  = mean(M1, dim);
    mn2  = mean(M2, dim);   
end

% Remove the mean from the data
M1_cnt = bsxfun(@minus, M1, mn1);
M2_cnt = bsxfun(@minus, M2, mn2);

% Compute the numerator and denominator separately
num   = sum(M1_cnt.*M2_cnt, dim);
denom = sqrt(sum(M1_cnt.^2, dim).*sum(M2_cnt.^2, dim));
    
R = num./denom;
    



