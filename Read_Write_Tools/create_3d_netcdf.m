function [] = create_3d_netcdf(fnme, dtainfo, varnme, varlong, varunits, times, lat, lon, time_units, chnks)   
    

    glbl_Atts = fieldnames(dtainfo);


    if nargin < 7
    	chnks = []; 
    end
    	
    ntimes = length(times);

    
    if size(lat, 2) > 1
        lat = lat(:, 1);
    end
    
    if size(lon, 2) > 1
        lon = lon(1, :)';
    end
    
    nlat   = length(lat);
    nlon   = length(lon);
        
    % Get the ID of the global attributes
    glob_id = netcdf.getConstant('NC_GLOBAL');

    % Set the parameters for NETCDF4-classic 
    cmode   = netcdf.getConstant('NETCDF4');
    cmode   = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));

    ncid        = netcdf.create(fnme, cmode);
    
    for i = 1:length(glbl_Atts)
        netcdf.putAtt(ncid, glob_id, glbl_Atts{i, :}, dtainfo.(glbl_Atts{i}));
    end

    
    time_dim_id = netcdf.defDim(ncid, 'time', ntimes);
    lat_dim_id  = netcdf.defDim(ncid, 'lat', nlat);
    lon_dim_id  = netcdf.defDim(ncid, 'lon', nlon);
    
    time_id     = netcdf.defVar(ncid, 'time', 'NC_FLOAT', time_dim_id);
    lat_id      = netcdf.defVar(ncid, 'lat', 'NC_FLOAT', lat_dim_id);
    lon_id      = netcdf.defVar(ncid, 'lon', 'NC_FLOAT', lon_dim_id);
    
    netcdf.putAtt(ncid, time_id, 'units', time_units);
    netcdf.putAtt(ncid, lat_id, 'units', 'degrees_north');
    netcdf.putAtt(ncid, lon_id, 'units', 'degrees_east');
    
    netcdf.putAtt(ncid, lat_id, 'standard_name', 'latitude');
    netcdf.putAtt(ncid, lon_id, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, time_id, 'standard_name', 'time');
    netcdf.putAtt(ncid, time_id, 'calendar', 'proleptic_gregorian');
    
    var_dims    = [lon_dim_id, lat_dim_id, time_dim_id];
    
    if iscell(varnme)
        for i = 1:length(varnme)
            var_id      = netcdf.defVar(ncid, varnme{i}, 'NC_FLOAT', var_dims);
            netcdf.defVarFill(ncid, var_id, false, NaN);
            
            if ~isempty(chnks)
                netcdf.defVarDeflate(ncid, var_id, false, true, 6);
                netcdf.defVarChunking(ncid,  var_id, 'CHUNKED', chnks);
            end
        end
    else
        var_id = netcdf.defVar(ncid, varnme, 'NC_FLOAT', var_dims);
        
        netcdf.putAtt(ncid, var_id, 'long_name', varlong);
        netcdf.putAtt(ncid, var_id, 'units', varunits);
        
        netcdf.defVarFill(ncid, var_id, false, NaN);
        if ~isempty(chnks)
            netcdf.defVarDeflate(ncid, var_id, false, true, 6);
            netcdf.defVarChunking(ncid,  var_id, 'CHUNKED', chnks);
        end
    end

    
    netcdf.endDef(ncid)

    ncwrite(fnme, 'time', single(times));
    ncwrite(fnme, 'lat', single(lat));
    ncwrite(fnme, 'lon', single(lon));      
end