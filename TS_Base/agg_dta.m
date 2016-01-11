function Out = agg_dta(fld, method)

if strcmp(method, 'mean')          
	Data_out = mean(fld, 1);
elseif strcmp(method, 'nanmean')
    Data_out = nanmean(fld, 1);
elseif strcmp(method, 'sum')
	Data_out = sum(fld, 1);
elseif strcmp(method, 'nansum')
    Data_out = nansum(fld, 1);
elseif strcmp(method, 'sum_squared')
    Data_out = sum(fld.^2, 1);
elseif strcmp(method, 'nansum_squared')
    Data_out = nansum(fld.^2, 1); 
elseif strcmp(method, 'rms')
    Data_out = rms(fld, 1);
elseif strcmp(method, 'nanrms')
    Data_out = nanrms(fld);
elseif strcmp(method, 'variance')
    Data_out = var(fld, 1);
elseif strcmp(method, 'nanvariance')
    Data_out = nanvar(fld, 1);
elseif strcmp(method, 'median')
    Data_out = median(fld, 1);
elseif strcmp(method, 'nanmedian')
    Data_out = nanmedian(fld, 1);
elseif strcmp(method, 'max')
    Data_out = max(fld, 1);
elseif strcmp(method, 'nanmax')
    Data_out = nanmax(fld, 1);
elseif strcmp(method, 'min')
    Data_out = min(fld, 1);
elseif strcmpt(method, 'nanmin')
    Data_out = nanmin(fld, 1);
end