function [mn, sig, g1, g2] = nanmoments(fld, dim);
%--------------------------------------------------------------------------
% The function computes the first four statistical moments from an input
% vector, matrix, or array. It also takes NaNs into account, as it removes
% them in advance. 
%--------------------------------------------------------------------------
% INPUT:
% - fld       Input vector, matrix, or array
% - dim       Dimension from which the moments should be derived 
%             (default: 1)
%--------------------------------------------------------------------------
% OUTPUT:
% - mn        Mean 
% - sig       Variance
% - g1        Skewness
% - g2        Curtosis
%--------------------------------------------------------------------------
% EXAMPLE:
% >> [mn, sig, g1, g2] = nanmoments(rand(100, 10, 10), 1);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Get the size of the field
sze_fld = size(fld);

% For counting the number of valid data, create a mask with ones.
mask             = ones(size(fld));
% Each element where isnan(fld) is set to 0
mask(isnan(fld)) = 0;
% Simply count the number of ones in the selected dimension
ndta             = sum(mask, dim);

% Compute the mean (first moment)
mn = nanmean(fld, dim);

% Remove the mean from the data
rep      = sze_fld(dim);
tmp      = ones(1, length(sze_fld));
tmp(dim) = rep;
    
fld_cnt  = fld - repmat(mn, tmp);


% Compute the population variance (second moment)
sig = 1./ndta.*nansum(fld_cnt.^2, dim);

% Compute the skewness (third moment)
g1 = (1./ndta.*nansum(fld_cnt.^3, dim))./...
                               ((1./ndta.*nansum(fld_cnt.^2, dim)).^(3/2));

% Compute the curtosis (fourth moment)
g2 = (1./ndta.*nansum(fld_cnt.^4, dim))./...
                               ((1./ndta.*nansum(fld_cnt.^2, dim)).^2) - 3;
                               
mn  = squeeze(mn);
sig = squeeze(sig);
g1  = squeeze(g1);
g2  = squeeze(g2);
