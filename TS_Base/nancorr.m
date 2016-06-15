function R = nancorr(fld1, fld2, dim)

if nargin < 3, dim = 1; end


fld1(isnan(fld2)) = NaN;
fld2(isnan(fld1)) = NaN;

% Compute the mean from both datasets
mn_fld1 = nanmean(fld1, dim);
mn_fld2 = nanmean(fld2, dim);
            
 % Remove the mean from the data
fld1_cnt = bsxfun(@minus, fld1, mn_fld1);
fld2_cnt = bsxfun(@minus, fld2, mn_fld2);

% Compute the sum of the products of the anomalies
num = nansum((fld1_cnt.*fld2_cnt), dim);

% Compute the products of the standard deviations
den = sqrt(nansum(fld1_cnt.^2, dim)).*sqrt(nansum(fld2_cnt.^2, dim));
           
R = num./den;