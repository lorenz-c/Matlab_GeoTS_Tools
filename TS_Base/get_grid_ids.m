function [ids dist] = get_grid_ids(lon, lat, lon_g, lat_g)


for i = 1:length(lon)
    lon_df = abs(lon_g - lon(i));
    lat_df = abs(lat_g - lat(i));
    
    [lat_dist(i, 1), lat_id(i, 1)] = min(lat_df);
    [lon_dist(i, 1), lon_id(i, 1)] = min(lon_df);
    
    metr_dist(i, 1) = haversine(lat(i), lon(i), lat_g(lat_id(i)), lon_g(lon_id(i)), 6378137, false);
end

ids  = [lon_id lat_id];
dist = [lon_dist, lat_dist metr_dist];


    