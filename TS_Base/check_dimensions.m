function [dta, nr_dims, out_dims, out_dim_nms] = check_dimensions(inpt, varnme, check_lat)
  

if nargin < 3, check_lat = false; end

% Read the number of dimensions of the actual variable
nr_dims = length(inpt.Variables.(varnme).dimensions);

% Get the names of the dimensions
for i = 1:nr_dims
    dim_nme{i} = inpt.Variables.(varnme).dimensions{i};
end

% Check for the time-dimension
if nr_dims == 3 || nr_dims == 4  
    try 
        time_dim_indx = find(ismember(dim_nme, 'time'));
    catch
        warning('Field is 3D but could not find the time-dimension')
    end
end

% Check for lat- and lon-dimensions
if nr_dims == 2 || nr_dims == 3 || nr_dims == 4
    try 
        lat_dim_indx = find(ismember(dim_nme, 'lat'));
    catch
        try 
            lat_dim_indx = find(ismember(dim_nme, 'latitude'));     
        catch
            warning('Could not find the lat-dimension')
        end
    end
    
    try 
        lon_dim_indx = find(ismember(dim_nme, 'lon'));
    catch
        try 
            lon_dim_indx = find(ismember(dim_nme, 'longitude'));    
        catch
            warning('Could not find the lon-dimension')
        end
    end
end

if nr_dims == 4
    dim_nme_tmp = dim_nme; 
    dim_nme_tmp([time_dim_indx, lat_dim_indx, lon_dim_indx]) = [];

    fourth_dim_indx = find(ismember(dim_nme, dim_nme_tmp{1}));
    
    out_dims    = [time_dim_indx, fourth_dim_indx, lat_dim_indx, lon_dim_indx];
    out_dim_nms = dim_nme(out_dims);
    
    dta = permute(inpt.Data.(varnme), out_dims);
     
elseif nr_dims == 3
    out_dims    = [time_dim_indx, lat_dim_indx, lon_dim_indx];
    out_dim_nms = dim_nme(out_dims);
    
    dta = permute(inpt.Data.(varnme), out_dims);
    
elseif nr_dims == 2
    out_dims = [lat_dim_indx, lon_dim_indx];
    out_dim_nms = dim_nme(out_dims);
    
    dta = permute(inpt.Data.(varnme), out_dims);
end

if check_lat == true
    if isfield(inpt.Data, 'lat')
        lat = inpt.Data.lat;
    elseif isfield(inpt.Data, 'latitude')
        lat = inpt.Data.latitude;
    else
        warning('Could not find latitude values')
    end
    
    if lat(1) < lat(end)
        
        lat = lat(end:-1:1);
        
        if nr_dims == 2
            dta = dta(end:-1:1, :);
        elseif nr_dims == 3
            dta = dta(:, end:-1:1, :);
        elseif nr_dims == 4
            dta = dta(:, :, end:-1:1, :);
        end
    end
end
    
    

    
    
    