function [Res, mnth_mn] = remsc_ts(inpt, tscale, remmn, frmt, clms)
% The function first estimates the seasonal cycle (i.e. the mean of
% January, February, ...) and removes it from the time-series in fld.
%--------------------------------------------------------------------------
% Input:        fld       [m1 x n]  Matrix which contains the input 
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
% References: B. D. Santer et al. (2000): Statistical significance of
% trends and trend differences in layer-average atmospheric temperature
% time series, Journal of Geophysical Research, 105: 7337-7356, 
% doi: 10.1029/1999JD901105
%--------------------------------------------------------------------------


if nargin < 4, frmt = 1; end
if nargin < 3, remmn = 0; end
if nargin < 2, tscale = 'monthly'; end

if strcmp(tscale, 'monthly')
	mnth_clm = 1;
	dta_clm  = 4;
elseif strcmp(tscale, 'daily')
    mnth_clm = 2;
    dta_clm  = 5;
elseif strcmp(tscale, 'none')
    dta_clm  = 1;
end

% The output will have the same dimensions as the input
Res = inpt;
   
% The loop goes over the twelve months

if frmt == 1
    nr_tstps = size(inpt, 1) - 1;
    
    if remmn == 1
        mn   = nanmean(inpt, 1);
        inpt = inpt - repmat(mn, nr_tstps, 1);
    end
    
    for i = 1:12
        
        % Find all months i in the input data
        mnth_indx = find(inpt(:, mnth_clm) == i);
    
        % Compute the monthly mean
        mnth_mn(i, :) = nanmean(inpt(mnth_indx, dta_clm:end));
        
        % Remove the seasonal cycle of month i from the corresponding columns
        Res(mnth_indx, dta_clm:end) = Res(mnth_indx, dta_clm:end) - ...
                                          ones(length(mnth_indx), 1)*mnth_mn(i, :);
    end
    
elseif frmt == 2
    nr_tstps = size(inpt.Time, 1);
    
    if remmn == 1
        inpt.Data = detrend(inpt.Data, 0);
    end
    
    for i = 1:12
        % Find all months i in the input data
        mnth_indx = find(inpt.Time(:, mnth_clm) == i);
        
        % Compute the monthly mean
        mnth_mn   = nanmean(inpt.Data(mnth_indx, :));
        
        % Remove the seasonal cycle of month i from the corresponding columns
        Res.Data(mnth_indx, :) = Res.Data(mnth_indx, :) - ...
                                        ones(length(mnth_indx), 1)*mnth_mn;
    end
    
elseif frmt == 3
    nr_tstps = size(inpt, 1);
    nr_years = nr_tstps/12;
    
	if remmn == 1
        mn   = nanmean(inpt, 1);
        inpt = inpt - repmat(mn, nr_tstps, 1);
	end

    
    
    for i = 1:12
        mnth_mn(i, :) = nanmean(inpt(i:12:end, :));
    end
    
    Res = Res - repmat(mnth_mn, nr_years, 1);
    
end
                                    
                                    
        
        
    
        
