function R = acf_mtrx(fld, l_max)

% Number of time-steps
nts = size(fld, 1);

% Compute the mean 
mn_fld  = ones(nts, 1)*nanmean(fld);

% Remove the mean from the field
fld_cnt = fld - mn_fld;

% Compute the unscaled co-variance
if l_max == 1
    sig_co = nansum(fld_cnt(2:end, :).*fld_cnt(1:end-1, :));
end 

% Compute the unscaled variance
sig_va = nansum(fld_cnt.*fld_cnt);

% Compute the autocorrelation
R = sig_co./sig_va;


