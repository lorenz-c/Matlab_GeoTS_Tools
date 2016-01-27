wrf    = netcdf2datastruct('http://schnepfe:8080/thredds/dodsC/WRF_Test/hydrology_d02_1971-2000all_WSM5_RAINall_dailytimestep.grd_CLM4.8_7km_regular_grid_GER.nc', false);
regnie = netcdf2datastruct('http://schnepfe:8080/thredds/dodsC/WRF_Test/REGNIE_1971-2000.nc_CLM4.8_7km_regular_grid_GER.nc', false);

regnie = renamevar(regnie, {'latitude', 'longitude'}, {'lat', 'lon'});

x  = wrf.Data.time;
T  = datestr(datenum(floor(x/10000), floor(mod(x,10000)/100), floor(mod(x,100)) ));
T2 = datenum(T);
T3 = datevec(T2);

wrf.DataInfo.title    = 'WRF';
regnie.DataInfo.title = 'Regnie';

[errs, TS_out, errs_out, anom_errs_out, obs_sc, mdl_sc] = batch_validation(regnie2, wrf, 'TOT_PREC', valid_options)