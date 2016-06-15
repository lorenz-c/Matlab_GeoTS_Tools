function ts_out = seltime(ts_in, time, unit)


vars    = fieldnames(ts_in.Variables),
hastime = istimevar(ts_in, vars); 

stat_vars          = vars(hastime == 0);
vars(hastime ~= 1) = [];

ts_out.DataInfo   = ts_in.DataInfo;
ts_out.Dimensions = ts_in.Dimensions;

for i = 1:length(stat_vars)
    ts_out.Variables.(stat_vars{i}) = ts_in.Variables.(stat_vars{i});
    ts_out.Data.(stat_vars{i})      = ts_in.Data.(stat_vars{i});
end


if strcmp(unit, 'minute')
    ids = find(ts_in.Data.time(:, 5) == time);
elseif strcmp(unit, 'hour')
    ids = find(ts_in.Data.time(:, 4) == time);
elseif strcmp(unit, 'day')
    ids = find(ts_in.Data.time(:, 3) == time);
elseif strcmp(unit, 'month')
    ids = find(ts_in.Data.time(:, 2) == time);
elseif strcmp(unit, 'year')
    ids = find(ts_in.Data.time(:, 1) == time);
end

time_ids = getdimpos(ts_in, vars, 'time');

for i = 1:length(vars)
    ts_out.Variables.(vars{i}) = ts_in.Variables.(vars{i});

    if time_ids(i) == 1
        if length(ts_in.Variables.(vars{i}).dimensions) == 1
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(ids, :);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 2
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(ids, :);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 3
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(ids, :, :);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 4
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(ids, :, :, :);
        end
    elseif time_ids(i) == 2
        if length(ts_in.Variables.(vars{i}).dimensions) == 2
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, ids);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 3
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, ids, :);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 4
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, ids, :, :);
        end
    elseif time_ids(i) == 3
        if length(ts_in.Variables.(vars{i}).dimensions) == 3
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, :, ids);
        elseif length(ts_in.Variables.(vars{i}).dimensions) == 4
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, :, ids, :);
        end
    elseif time_ids(i) == 4
        if length(ts_in.Variables.(vars{i}).dimensions) == 4
            ts_out.Data.(vars{i}) = ts_in.Data.(vars{i})(:, :, :, ids);
        end
    end
end

ts_out.TimeStamp = ts_in.TimeStamp(ids, :);


if ~isinf(ts_out.Dimensions.time)
    ts_out.Dimensions.time = length(ts_out.TimeStamp);
end
