function [ts_out varargout] = ts_average(ts_in, tres_out, method, vars)
% Compute temporal averages of a time-series. The temporal resolution of
% the output time-series is defined by the tres_out parameter. Note that
% the function does not perform an interpolation, i.e. monthly data can NOT
% be downscaled to e.g. daily data! Furthermore, the function accepts
% different "methods" for the aggregation. If not all variables in an input
% datastructure should be aggregated, specific variables can be selected by
% the "vars"-parameter. 
% The function always assumes the first dimension to be the
% "time"-dimension, i.e. the aggregation is always performed over time.
% Depending on the selected tres_out, the output will contain either a
% "climatology_bounds" or "time_bounds" variable and the corresponding
% attributes. Thus, the output is conform with the CF-conventions which
% makes it easier to understand and interpret the data. More information 
% about the climatology- and time-bounds meta-data as well as the 
% cell_methods attributes can be found at http://cfconventions.org. Some 
% examples are presented in chapter 7.3 of the version 1.7 document.
% -------------------------------------------------------------------------
% Input:        
% - ts_in       Input datastructure
% - tres_out    String which defines the time scale of the output; valid 
%               options are: annual_lt (long-term mean), seasonal_lt 
%               (long-term MAM, JJA, SON, DJF mean), monthly_lt (long-term 
%               monthly mean -> mean annual cycle), annual (annual mean), 
%               seasonal (MAM, JJA, SON, DJF mean for every year -> 3-month 
%               block average), monthly (monthly mean from e.g. daily 
%               data), daily (daily mean from e.g. hourly data). 
% - method      String which defines the aggregation method; valid options
%               are: mean, nanmean (NaNs are neglected; default), sum, 
%               nansum, sum_squared, nansum_squared, rms, nanrms, variance,
%               nanvariance, median, nanmedian, max, nanmax, min, nanmin
% - vars        Can be set to "all" (default) -> aggregation for all 
%               variables with a "time"-dimension. If only selected 
%               variables should be aggregated, they should be provided in 
%               a cell array (i.e. vars = {'prec', 'temp', ...}.
% -------------------------------------------------------------------------
% Output:       
% - otpt      	Output datastructure
% -------------------------------------------------------------------------
% Author: Christof Lorenz
% Date:   July 2011
% -------------------------------------------------------------------------
% Uses: nanrms.m
% -------------------------------------------------------------------------
% Updates: - 01.04.2015: Brush up code, added some comments (C. Lorenz)
%          - 04.12.2015: Added support for CF-datastructures; added more 
%                        methods including min, max, median,...; added 
%                        comments, etc. (C. Lorenz)
% -------------------------------------------------------------------------

% Check input variables and set default parameter
if nargin < 4 | isempty(vars), vars = 'all'; end
if nargin < 3 | isempty(method), method = 'nanmean'; end

% Use the inpt datastructure as basis for the output
if isfield(ts_in, 'DataInfo')
    ts_out.DataInfo   = ts_in.DataInfo;
end
ts_out.Dimensions = ts_in.Dimensions;

if strcmp(vars, 'all')
    % Get the variables in the file
    vars = fieldnames(ts_in.Variables);
    
    % Check for fixed variables
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
        
        % Copy the fixed variables to the output (but don't do that for the
        % "time"-data as the time-vector of the output will look different
        if isfixed(i) == 1
            ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
            
            if ~strcmp(vars{i}, 'time')
                ts_out.Data.(vars{i}) = ts_in.Data.(vars{i});
            end
        end    
    end
    
    % Remove all fixed variables
    vars(isfixed == 1) = [];
    
    % Check for variables with the time-dimension
    for i = 1:length(vars)
        if isfield(ts_in.Variables.(vars{i}), 'dimensions')
            if find(ismember(ts_in.Variables.(vars{i}).dimensions,'time'))
                hastime(i) = 1;
            else
                hastime(i) = 0;
            end
        else
            hastime(i) = 0;
        end
        
        % Copy the non-fixed static variables to the output
        if hastime(i) == 0
            ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
            if isfield(ts_out.Data, vars{i})
                ts_out.Data.(vars{i})  = ts_in.Data.(vars{i});
            end
        end  
    end
    % For the following computations remove all variables which do not have 
    % a "time"-dimension
    vars(hastime == 0) = [];    
    
else
    
    % If only specific variables should be written to the output, we still 
    % need the "fixed" variables (lat, lon, ...)
    vars_old = fieldnames(ts_in.Variables);
    
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars_old{i});
        
        if isfixed(i) == 1
            ts_out.Variables.(vars_old{i}) = ts_in.Variables.(vars_old{i});
            
            if ~strcmp(vars{i}, 'time')
                ts_out.Data.(vars_old{i})  = ts_in.Data.(vars_old{i});
            end
        end    
    end   
end


% -------------------------------------------------------------------------
%                       Long term average (annual_lt)
% -------------------------------------------------------------------------
if strcmp(tres_out, 'annual_lt') 
        
    for i = 1:length(vars)
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        if tme_indx == 1
            ts_out.Data.(vars{i}) = agg_ts(ts_in.Data.(vars{i}), method);
            nr_nan.(vars{i})      = agg_ts(nan_mask, 'sum');
        elseif tme_indx == 2
            ts_out.Data.(vars{i}) = agg_ts(ts_in.Data.(vars{i}), method, 2);
            nr_nan.(vars{i})      = agg_ts(nan_mask, 'sum', 2);
        end
        
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
    
    end
    
    ts_out.Dimensions.time            = 1;
    ts_out.Dimensions.nv              = 2;

    ts_out.Variables.time.climatology = 'climatology_bounds';
    ts_out.Data.time                  = ts_in.Data.time(1, :);
    ts_out.TimeStamp                  = datenum(ts_out.Data.time);
    
    ts_out.Variables.climatology_bounds.dimensions = {'time', 'nv'};
    ts_out.Variables.climatology_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
    ts_out.Data.climatology_bounds                 = ...
                           [ts_in.Data.time(1, :) ts_in.Data.time(end, :)];
                                      
                       
% -------------------------------------------------------------------------
%                   Seasonal long term average (seasonal_lt)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'seasonal_lt')
    
    for i = 1:length(vars)
 
        sze_dta = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = getdimpos(ts_in, vars{i}, 'time');
                                                                 
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        
        % Create a vector which contains the months where data is available
        mnths = ts_in.Data.time(:, 2);
    
        % Loop over the four seasons
        for j = 1:4
            
            nan_mask = zeros(size(ts_in.Data.(vars{i})));
            nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
            if j == 4
                indx_1 = find(mnths == 12);
                indx_2 = find(mnths == 1);
                indx_3 = find(mnths == 2);            
            else
                indx_1 = find(mnths == (j+1)*3-3);
                indx_2 = find(mnths == (j+1)*3-2);
                indx_3 = find(mnths == (j+1)*3-1);
            end
            
        
            indx = [indx_1; indx_2; indx_3];
            
            bnds_first(j, :) = ts_in.Data.time(indx(1), :);
            bnds_last(j, :)  = ts_in.Data.time(indx(end), :);
            
            if j == 1
                yr_first = ts_in.Data.time(indx(1), 1);
            end
            
            if length(sze_dta) == 2
                if tme_indx == 1
                    ts_out.Data.(vars{i})(j, :) = ...
                             agg_ts(ts_in.Data.(vars{i})(indx, :), method);
                    nr_nan.(vars{i})(j, :) = ...
                                          agg_ts(nan_mask(indx, :), 'sum');
                elseif tme_indx == 2
                    ts_out.Data.(vars{i})(:, j) = ...
                          agg_ts(ts_in.Data.(vars{i})(:, indx), method, 2);
                    nr_nan.(vars{i})(:, j) = ...
                                       agg_ts(nan_mask(:, indx), 'sum', 2);
                end
            elseif length(sze_dta) == 3
                if tme_indx == 1
                    ts_out.Data.(vars{i})(j, :, :) = ...
                          agg_ts(ts_in.Data.(vars{i})(indx, :, :), method);
                    nr_nan.(vars{i})(j, :, :) = ...
                                       agg_ts(nan_mask(indx, :, :), 'sum');
                else
                    error('Wrong dimensions')
                end
            end
                  
            clear indx_1 indx_2 indx_3
        end

        % Add "climatology" meta-data to the output
        ts_out.Dimensions.time  = 4;
        ts_out.Dimensions.nv = 2;

        ts_out.Variables.time.climatology = 'climatology_bounds';
        ts_out.Data.time                  = ...
                          [yr_first 4  15 0 0 0; yr_first 7  15 0 0 0;  ...
                           yr_first 10 15 0 0 0; yr_first + 1 01 15 0 0 0];
                                         
        ts_out.TimeStamp  = datenum(ts_out.Data.time);
    
        ts_out.Variables.climatology_bounds.dimensions = {'time', 'nv'};
        ts_out.Variables.climatology_bounds.units = 'yyyy-mm-dd HH:MM:SS';
        ts_out.Data.climatology_bounds = [bnds_first bnds_last];
    end
   
    
% -------------------------------------------------------------------------
%                    Monthly long term average (monthly_lt)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'monthly_lt')
    
    % Create a vector which contains the months where data is available
    mnths = ts_in.Data.time(:, 2);
    
    for i = 1:12
        % Search for all months in the data which correspond to the ijth 
        % month:
        indx{i} = find(mnths == i);
 
        indx_first(i, 1) = indx{i}(1);
        indx_last(i, 1)  = indx{i}(end);
        
        date_first(i, :) = ts_in.Data.time(indx_first(i), :);
        date_last(i, :)  = ts_in.Data.time(indx_last(i), :);    
    end
    
    for i = 1:length(vars)

        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
                                                                 
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                         ['time: ', method, ' over years'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
        % Loop over twelve months 
        for j = 1:12          
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp  = ts_in.Data.(vars{i})(indx{j}, :);
                    tmp2 = nan_mask(indx{j}, :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp  = ts_in.Data.(vars{i})(:, indx{j});
                    tmp2 = nan_mask(:, indx{j});
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                if tme_indx == 1
                    tmp  = ts_in.Data.(vars{i})(indx{j}, :, :);
                    tmp2 = nan_mask(indx{j}, :, :);
                    ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
                else
                    error('Wrong dimensions!')
                end
            end           
        end
        
        ts_out.Dimensions.time            = 12;
        ts_out.Dimensions.nv              = 2;

        ts_out.Variables.time.climatology = 'climatology_bounds';
        ts_out.Data.time                  = date_first;
                                         
        ts_out.TimeStamp                  = datenum(ts_out.Data.time);
    
        ts_out.Variables.climatology_bounds.dimensions = {'time', 'nv'};
        ts_out.Variables.climatology_bounds.units = 'yyyy-mm-dd HH:MM:SS';
        ts_out.Data.climatology_bounds = [date_first date_last];
        
    end
    
elseif strcmp(tres_out, 'daily_ltm')
    
    dys = 1:365;
    
    for i = 1:length(ts_in.Data.time)
        if ts_in.Data.time(i, 2) == 2 && ts_in.Data.time(i, 3) == 29
            doy(i) = 59;
        else
            frstjan = datetime([ts_in.Data.time(i, 1), 1,1]);
            actdte  = datetime(ts_in.Data.time(i, :));
            doy(i)  = days(actdte - frstjan) + 1;
        end
    end
 
            
    for i = 1:365
        
        indx{i} = find(doy == i);
 
        indx_first(i, 1) = indx{i}(1);
        indx_last(i, 1)  = indx{i}(end);
        
        date_first(i, :) = ts_in.Data.time(indx_first(i), :);
        date_last(i, :)  = ts_in.Data.time(indx_last(i), :);    
    end
    
    for i = 1:length(vars)

        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
                                                                 
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                         ['time: ', method, ' over years'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
        % Loop over twelve months 
        for j = 1:365         
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp  = ts_in.Data.(vars{i})(indx{j}, :);
                    tmp2 = nan_mask(indx{j}, :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp  = ts_in.Data.(vars{i})(:, indx{j});
                    tmp2 = nan_mask(:, indx{j});
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                if tme_indx == 1
                    tmp  = ts_in.Data.(vars{i})(indx{j}, :, :);
                    tmp2 = nan_mask(indx{j}, :, :);
                    ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
                else
                    error('Wrong dimensions!')
                end
            end           
        end
        
        ts_out.Dimensions.time            = 12;
        ts_out.Dimensions.nv              = 2;

        ts_out.Variables.time.climatology = 'climatology_bounds';
        ts_out.Data.time                  = date_first;
                                         
        ts_out.TimeStamp                  = datenum(ts_out.Data.time);
    
        ts_out.Variables.climatology_bounds.dimensions = {'time', 'nv'};
        ts_out.Variables.climatology_bounds.units = 'yyyy-mm-dd HH:MM:SS';
        ts_out.Data.climatology_bounds = [date_first date_last];
        
    end
    
    
% -------------------------------------------------------------------------
%                          Annual average (annual)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'annual')
    
	yrs = unique(ts_in.Data.time(:, 1));
    
    for i = 1:length(yrs)
        indx = find(ts_in.Data.time(:, 1) == yrs(i));
        
        indx_first(i, 1) = indx(1);
        indx_last(i, 1)  = indx(end);
        
        date_first(i, :) = ts_in.Data.time(indx(1), :);
        date_last(i, :)  = ts_in.Data.time(indx(end), :); 
    end
    
    for i = 1:length(vars)
        
        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        % Change the size of the time-dimension in the output                                                        
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                 ['time: ', method, ' (interval: 1 year)'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
        for j = 1:length(indx_first)
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp = ...
                       ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp = ...
                       ts_in.Data.(vars{i})(:, indx_first(j):indx_last(j));
                    tmp2 = nan_mask(:, indx_first(j):indx_last(j));
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                if tme_indx == 1
                    tmp = ...
                    ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :, :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :, :);
                    ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
                else
                    error('Wrong dimensions!')
                end
            end              
        end   
    end
    
    ts_out.Dimensions.time       = Inf;
    ts_out.Dimensions.nv         = 2;

    ts_out.Variables.time.bounds = 'time_bounds';
    ts_out.Data.time             = [date_first(:, 1), ...
                                    ones(length(yrs), 2), ...
                                    zeros(length(yrs), 3)];
                                         
    ts_out.TimeStamp             = datenum(ts_out.Data.time);
    
	ts_out.Variables.time_bounds.dimensions = {'time', 'nv'};
	ts_out.Variables.time_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
	ts_out.Data.time_bounds                 = [date_first date_last];
        
        
% -------------------------------------------------------------------------
%                        Seasonal average (seasonal)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'seasonal')
    
    % Get the first and last year of the time-series
    syr = ts_in.Data.time(1, 1);
    eyr = ts_in.Data.time(end, 1);
    
    % Vector of unique years
    yrs_unique = unique(ts_in.Data.time(:, 1));
    yrs        = ts_in.Data.time(:, 1);
    mnths      = ts_in.Data.time(:, 2);
    
    k = 1;
    for i = 1:length(yrs_unique)
        for j = 1:4
            if j == 1 && i > 1
                indx_1 = find(mnths == 12 & yrs == yrs_unique(i)-1);
                indx_2 = find(mnths == 2  & yrs == yrs_unique(i));  
            elseif j ~= 1
                indx_1 = find(mnths == j*3-3 & yrs == yrs_unique(i));
                indx_2 = find(mnths == j*3-1 & yrs == yrs_unique(i));
            else
                indx_1 = [];
                indx_2 = [];     
            end

            if ~isempty(indx_1) & ~isempty(indx_2)
                indx_first(k, 1) = indx_1(1);
                indx_last(k, 1)  = indx_2(end);
                if j == 1 && i > 1
                    ts_out.Data.time(k, :) = [yrs_unique(i) ...
                                                               1 15 0 0 0];
                elseif j ~= 1
                    ts_out.Data.time(k, :) = [yrs_unique(i) ...
                                                     ((j-1)*3+1) 15 0 0 0];
                end
                k = k + 1;
            end 
        end
    end

    date_first = ts_in.Data.time(indx_first, :);
    date_last  = ts_in.Data.time(indx_last, :); 
            
    for i = 1:length(vars)
        
        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        % Change the size of the time-dimension in the output                                                        
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                               ['time: ', method, ' (interval: 1 season)'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
        for j = 1:length(indx_first)
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp = ...
                       ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp = ...
                       ts_in.Data.(vars{i})(:, indx_first(j):indx_last(j));
                    tmp2 = nan_mask(:, indx_first(j):indx_last(j));
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                if tme_indx == 1
                    tmp = ...
                    ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :, :);
                    ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp, 'sum');
                else
                    error('Wrong dimensions!')
                end       
            end              
        end   
    end
    
    ts_out.Dimensions.time       = Inf;
    ts_out.Dimensions.nv         = 2;

	ts_out.Variables.time.bounds = 'time_bounds';
                                   
	ts_out.TimeStamp              = datenum(ts_out.Data.time);
    
    ts_out.Variables.time_bounds.dimensions = {'time', 'nv'};
    ts_out.Variables.time_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
    ts_out.Data.time_bounds                 = [date_first date_last];
        
    ts_out.TimeStamp = datenum(ts_out.Data.time);
    
    
% -------------------------------------------------------------------------
%                         Monthly average (monthly)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'monthly')
    
    % Vector of unique years
    yrs = unique(ts_in.Data.time(:, 1));
    
    k = 1;
    for i = 1:length(yrs)
        for j = 1:12
            indx = find(ts_in.Data.time(:, 1) == yrs(i) & ...
                        ts_in.Data.time(:, 2) == j);
                    
            if ~isempty(indx)
                indx_first(k, 1) = indx(1);
                indx_last(k, 1)  = indx(end);
                
                ts_out.Data.time(k, :) = [yrs(i) j 15 0 0 0];
                
                k = k + 1;
            end
        end
    end
    
    date_first = ts_in.Data.time(indx_first, :);
    date_last  = ts_in.Data.time(indx_last, :); 
    
    for i = 1:length(vars)
       
        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        % Change the size of the time-dimension in the output                                                        
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                ['time: ', method, ' (interval: 1 month)'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;
        
        for j = 1:length(indx_first)
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp = ...
                       ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp = ...
                       ts_in.Data.(vars{i})(:, indx_first(j):indx_last(j));
                    tmp2 = nan_mask(:, indx_first(j):indx_last(j));
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                tmp = ...
                    ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :, :);
                tmp2 = nan_mask(indx_first(j):indx_last(j), :, :);
                ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
            end              
        end   
    end
    
    ts_out.Dimensions.time       = Inf;
    ts_out.Dimensions.nv         = 2;

	ts_out.Variables.time.bounds = 'time_bounds';
                                   
	ts_out.TimeStamp              = datenum(ts_out.Data.time);
    
    ts_out.Variables.time_bounds.dimensions = {'time', 'nv'};
    ts_out.Variables.time_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
    ts_out.Data.time_bounds                 = [date_first date_last];
        
    ts_out.TimeStamp = datenum(ts_out.Data.time);
    
% -------------------------------------------------------------------------
%                        Weekly (daily)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'weekly')
    
    % Vector of unique years
    yrs = unique(ts_in.Data.time(:, 1));
    
    l = 1;
    for i = 1:length(yrs)
        
        yr_indx   = find(ts_in.Data.time(:, 1) == yrs(i));
        t_vec_act = ts_in.TimeStamp(yr_indx);
        
        % Get the weeks of the year
        weeks     = weeknum(t_vec_act);
        weeks_uni = unique(weeks);
        
        for j = 1:length(weeks_uni)
            
            week_indx    = find(weeks == weeks_uni(j));
            sdte         = t_vec_act(week_indx(1));
            edte         = t_vec_act(week_indx(end));
            
            indx_first(l, 1) = find(ts_in.TimeStamp == sdte);
            indx_last(l, 1)  = find(ts_in.TimeStamp == edte);
            
            ts_out.Data.time(l, :) = datevec(sdte);
            
            l = l + 1;
        end 
    end
    
    date_first = ts_in.Data.time(indx_first, :);
    date_last  = ts_in.Data.time(indx_last, :); 
    
    for i = 1:length(vars)  
        % Get the size of the input data
        sze_dta  = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        % Change the size of the time-dimension in the output                                                        
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                ['time: ', method, ' (interval: 7 days)'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;

        for j = 1:length(indx_first)
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp = ...
                       ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp = ...
                       ts_in.Data.(vars{i})(:, indx_first(j):indx_last(j));
                    tmp2 = nan_mask(:, indx_first(j):indx_last(j));
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                tmp = ...
                    ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :, :);
                tmp2 = nan_mask(indx_first(j):indx_last(j), :, :);
                ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
            end              
        end   
    end
    
    ts_out.Dimensions.time       = Inf;
    ts_out.Dimensions.nv         = 2;

	ts_out.Variables.time.bounds = 'time_bounds';
                                   
	ts_out.TimeStamp              = datenum(ts_out.Data.time);
    
    ts_out.Variables.time_bounds.dimensions = {'time', 'nv'};
    ts_out.Variables.time_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
    ts_out.Data.time_bounds                 = [date_first date_last];
        
    ts_out.TimeStamp = datenum(ts_out.Data.time);
    
    
% -------------------------------------------------------------------------
%                        Daily average (daily)
% -------------------------------------------------------------------------
elseif strcmp(tres_out, 'daily')
    
    % Vector of unique years
    yrs = unique(ts_in.Data.time(:, 1));
    
    l = 1;
    for i = 1:length(yrs)
        for j = 1:12
            nrd = eomday(yrs(i), j);
            for k = 1:nrd
                indx = find(ts_in.Data.time(:, 1) == yrs(i) & ...
                            ts_in.Data.time(:, 2) == j & ...
                            ts_in.Data.time(:, 3) == k);
                    
                if ~isempty(indx)
                    indx_first(l, 1) = indx(1);
                    indx_last(l, 1)  = indx(end);
                
                    ts_out.Data.time(l, :) = [yrs(i) j k 0 0 0];
                
                    l = l + 1;
                end
            end
        end
    end
    
    date_first = ts_in.Data.time(indx_first, :);
    date_last  = ts_in.Data.time(indx_last, :); 
    
    for i = 1:length(vars)  
        % Get the size of the input data
        sze_dta    = size(ts_in.Data.(vars{i}));
        
        % Get the "position" of the time-dimension
        tme_indx = ...
            find(ismember(ts_in.Variables.(vars{i}).dimensions, 'time') ...
                                                                     == 1);
        % Change the size of the time-dimension in the output                                                        
        sze_dta(tme_indx) = length(indx_first);
        
        % Copy the variable's meta-data to the output variable
        ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});
        % Add the cell-methods attribute 
        ts_out.Variables.(vars{i}).cell_methods = ...
                                ['time: ', method, ' (interval: 1 day)'];
        % Create an empty array for the data                  
        ts_out.Data.(vars{i}) = NaN(sze_dta);
        % Create an empty array for counting the number of NaNs      
        nan_mask = zeros(size(ts_in.Data.(vars{i})));
        % nan_mask == 1 where ts_in.Data.(vars{i}) has missing values.
        nan_mask(isnan(ts_in.Data.(vars{i}))) = 1;

        for j = 1:length(indx_first)
            if length(sze_dta) == 2
                if tme_indx == 1
                    tmp = ...
                       ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :);
                    tmp2 = nan_mask(indx_first(j):indx_last(j), :);
                    ts_out.Data.(vars{i})(j, :) = agg_ts(tmp, method);
                    nr_nan.(vars{i})(j, :)      = agg_ts(tmp2, 'sum');
                elseif tme_indx == 2
                    tmp = ...
                       ts_in.Data.(vars{i})(:, indx_first(j):indx_last(j));
                    tmp2 = nan_mask(:, indx_first(j):indx_last(j));
                    ts_out.Data.(vars{i})(:, j) = agg_ts(tmp, method, 2);
                    nr_nan.(vars{i})(:, j)      = agg_ts(tmp2, 'sum', 2);
                end
            elseif length(sze_dta) == 3
                tmp = ...
                    ts_in.Data.(vars{i})(indx_first(j):indx_last(j), :, :);
                tmp2 = nan_mask(indx_first(j):indx_last(j), :, :);
                ts_out.Data.(vars{i})(j, :, :) = agg_ts(tmp, method);
                nr_nan.(vars{i})(j, :, :)      = agg_ts(tmp2, 'sum');
            end              
        end   
    end
    
    ts_out.Dimensions.time       = Inf;
    ts_out.Dimensions.nv         = 2;

	ts_out.Variables.time.bounds = 'time_bounds';
                                   
	ts_out.TimeStamp              = datenum(ts_out.Data.time);
    
    ts_out.Variables.time_bounds.dimensions = {'time', 'nv'};
    ts_out.Variables.time_bounds.units      = 'yyyy-mm-dd HH:MM:SS';
    ts_out.Data.time_bounds                 = [date_first date_last];
        
    ts_out.TimeStamp = datenum(ts_out.Data.time);
end

% Put nr_nan to varargout
% varargout{1} = nr_nan;


% Update the file history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
    '; MATLAB TS-Tools: ts_average.m -> ', method, ' average'];    
if isfield(ts_out, 'DataInfo')        
    ts_out.DataInfo.history = ...
                        sprintf([new_hist, ' \n', ts_in.DataInfo.history]);
else
    ts_out.DataInfo.history = sprintf([new_hist, ' \n']);
end

end
            
      
function Data_out = agg_ts(fld, method, dim)
    if nargin < 3
        dim = 1;
    end
    
    if strcmp(method, 'mean')          
        Data_out = mean(fld, dim);
    elseif strcmp(method, 'nanmean')
        Data_out = nanmean(fld, dim);
    elseif strcmp(method, 'sum')
        Data_out = sum(fld, dim);
    elseif strcmp(method, 'nansum')
        Data_out = nansum(fld, dim);
    elseif strcmp(method, 'sum_squared')
        Data_out = sum(fld.^2, dim);
    elseif strcmp(method, 'nansum_squared')
        Data_out = nansum(fld.^2, dim); 
    elseif strcmp(method, 'rms')
        Data_out = rms(fld, dim);
    elseif strcmp(method, 'nanrms')
        Data_out = nanrms(fld);
    elseif strcmp(method, 'std')
        Data_out = std(fld, dim);
    elseif strcmp(method, 'nanstd')
        Data_out = nanstd(fld, dim);
    elseif strcmp(method, 'variance')
        Data_out = var(fld, dim);
    elseif strcmp(method, 'nanvariance')
        Data_out = nanvar(fld, dim);
    elseif strcmp(method, 'median')
        Data_out = median(fld, dim);
    elseif strcmp(method, 'nanmedian')
        Data_out = nanmedian(fld, dim);
    elseif strcmp(method, 'max')
        Data_out = max(fld, dim);
    elseif strcmp(method, 'nanmax')
        Data_out = nanmax(fld, dim);
    elseif strcmp(method, 'min')
        Data_out = min(fld, dim);
    elseif strcmpt(method, 'nanmin')
        Data_out = nanmin(fld, dim);
    end
end


















            
        