function TS_out = find_sim_tstps(ref_series, TS_in, ref_time, window)

if nargin < 4, window = 1; end

% STILL SOME THINGS TO DO!!!!!!!!!!!!!!!!!!!!!

% Get all variables with a time-dimensions
vars              = fieldnames(TS_in.Variables);
istime            = istimevar(TS_in, vars);
vars(istime == 0) = [];

time_pos = getdimpos(TS_in, vars, 'time');


inpt_num = TS_in.TimeStamp;

TS_out = TS_in;

for i = 1:length(vars)
    sze_dta_old = size(TS_in.Data.(vars{i}));
    sze_dta_new = sze_dta_old;
    sze_dta_new(time_pos(i)) = size(ref_series, 1);
    
    TS_out.Data.(vars{i}) = NaN(sze_dta_new);
end


for i = 1:size(ref_series, 2)
    df = bsxfun(@minus, ref_series(:, i), inpt_num');
    df = abs(df);
    
    if window > 1
        [srted, idx] = sort(df, 2, 'ascend');
    else
        [srted, idx] = min(df, [], 2);
    end
    
    
    for j = 1:length(vars)
        if ~strcmp(vars{j}, 'time')
            if time_pos(j) == 1
                TS_tmp = TS_in.Data.(vars{j})(idx(:, 1), i);
                TS_tmp(isnan(srted(idx, 1))) = NaN;
                TS_out.Data.(vars{j})(:, i)  = TS_tmp;
            elseif time_pos(j) == 2
                TS_tmp = TS_in.Data.(vars{j})(i, idx(:, 1));
                TS_tmp(1, isnan(srted(:, 1))) = NaN;
                if window > 1
                    for k = 2:window
                        % Look for NaN-values in the output and get the
                        % second (third, fourth, ...) smallest value in the
                        % row
                        nan_indx         = find(isnan(TS_tmp));
                        TS_tmp(nan_indx) = TS_in.Data.(vars{j})(i, idx(nan_indx, k));
                        TS_tmp(isnan(TS_in.Data.(vars{j})(i, idx(nan_indx, k)))) = NaN;
                    end
                end
                TS_tmp(isnan(ref_series(:, i))) = NaN;
                TS_out.Data.(vars{j})(i, :) = TS_tmp;
            end
        end
    end

end


    
TS_out.Data.time = ref_time;
TS_out.TimeStamp = datenum(ref_time);


TS_out.Data.match_time = ref_series - datenum('1800-01-01 00:00:00');;

if time_pos(end) == 2
    TS_out.Data.match_time = TS_out.Data.match_time';
end

TS_out.Variables.match_time.long_name  = 'Exact time of matching time steps';
TS_out.Variables.match_time.units      = 'days since 1800-01-01 00:00:00';
TS_out.Variables.match_time.dimensions = TS_out.Variables.(vars{end}).dimensions;
TS_out.Variables.match_time.FillValue  = NaN;
