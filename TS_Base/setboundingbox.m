function otpt = setmask(inpt, bbox, var)

% Copy the global metadata to the output
otpt.Datainfo = inpt.DataInfo;

if size(inpt.Data.lat, 2) > 1
    if size(inpt.Data.lat, 1) ~= 1
        lat = inpt.Data.lat(:, 1);
    else
        lat = inpt.Data.lat';
    end
else
    lat = inpt.Data.lat;
end

if size(inpt.Data.lon, 2) > 1
    if size(inpt.Data.lon, 1) ~= 1
        lon = inpt.Data.lon(:, 1);
    else
        lon = inpt.Data.lon';
    end
else
    lon = inpt.Data.lon;
end

d_lon_left   = abs(inpt.Data.lon - bbox(1));
d_lon_right  = abs(inpt.Data.lon - bbox(2));
d_lat_top    = abs(inpt.Data.lat - bbox(3));
d_lat_bottom = abs(inpt.Data.lat - bbox(4));

[~, left_indx]   = min(d_lon_left);
[~, right_indx]  = min(d_lon_right);
[~, top_indx]    = min(d_lat_top);
[~, bottom_indx] = min(d_lat_bottom);





% Copy all fixed variables to the output
vars = fieldnames(inpt.Variables);
for i = 1:length(vars)
    isfixed(i) = isfixedvar(vars{i});
end

vars_fixed         = vars(isfixed == 1);
vars(isfixed == 1) = [];

otpt               = copyvars(otpt, inpt, vars_fixed);

% Check for gridded variables
isgrid = isgridvar(inpt, vars);

% Remove variables which do not have the lat- and lon-dimension
vars_nogrid = vars(isgrid == 0);

otpt = copyvars(otpt, inpt, vars_nogrid);

vars(isgrid == 0) = [];


for i = 1:length(vars)
    % Get the dimension
    dims = inpt.Variables.(vars{i}).dimensions;
    
    if length(dims) == 3
       % keyboard
        mask_3d                = reshape(mask, [1 size(mask)]);

        otpt.Data.(vars{i}) = bsxfun(@times, otpt.Data.(vars{i}), mask_3d);
    elseif length(dims) == 2
        otpt.Data.(vars{i}) = otpt.Data.(vars{i}).*mask;
    end
end

% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                           '; MATLAB TS-Tools: setmask.m'];
if isfield(inpt.DataInfo, 'history')           
    otpt.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                   inpt.DataInfo.history]);
else
    otpt.DataInfo.history = new_hist;
end

