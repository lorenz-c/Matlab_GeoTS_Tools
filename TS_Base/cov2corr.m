function C = cov2corr(Q)
%--------------------------------------------------------------------------
% The function transforms a covariance matrix Q in a correlation matrix C
% by dividing each matrix element by the correct product of the main
% diagonal elements (i.e. the variances and co-variances).
%--------------------------------------------------------------------------
% INPUT:
% - Q         Covariance matrix
%--------------------------------------------------------------------------
% OUTPUT:
% - C         Correlation matrix with the same dimensions as the input
%             matrix
%--------------------------------------------------------------------------
% EXAMPLE:
% >> C = cov2corr(Q);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Get the inverse square root of the diagonal elements -> inverse standard
% deviation
Vars = diag(diag(Q).^(-1/2));
% Compute the correlation matrix by multiplying the coviariance from both
% sides with the appropriate inverse standard deviations
C    = Vars*Q*Vars;

