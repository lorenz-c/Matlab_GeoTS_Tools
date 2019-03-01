function P = empirical_probab_forecast(inpt_dir, varnme, climatology, selected_percentiles, maxlead, fnme_out)
   
if nargin < 5
    maxlead = 6;
end

if nargin < 4
    selected_percentiles = [0 0.33 0.66 1]';
end

% Get the filenames of all files in the current directory
fles     = dir([inpt_dir, '/*', varnme, '*.nc']);


% Get the time of the forecast in the current directory
tme      = ncread([fles(1).folder, '/', fles(1).name], 'time');
tme_unit = ncreadatt([fles(1).folder, '/', fles(1).name], 'time', 'units');
tme_abs  = reldate2absdate(tme, tme_unit);

% Get the latitudes and longitudes
lat      = ncread([fles(1).folder, '/', fles(1).name], 'lat');
lon      = ncread([fles(1).folder, '/', fles(1).name], 'lon');

% Calculate the lenghts of the dimensions
nlat    = length(lat);
nlon    = length(lon);
ncategs = length(selected_percentiles) - 1;
nfles   = length(fles);

% Create empty matrices for all output variables
probab      = NaN(maxlead + 1, ncategs, nlon, nlat);
categ_mn    = NaN(maxlead + 1, ncategs, nlon, nlat);
categ_std   = NaN(maxlead + 1, ncategs, nlon, nlat);
flgs_max    = NaN(maxlead + 1,          nlon, nlat);
ens_mn      = NaN(maxlead + 1,          nlon, nlat);
ens_std     = NaN(maxlead + 1,          nlon, nlat);
ens_iqr     = NaN(maxlead + 1,          nlon, nlat);
ens_spread  = NaN(maxlead + 1,          nlon, nlat);
diff_climat = NaN(maxlead + 1,          nlon, nlat);
climat      = NaN(maxlead + 1,          nlon, nlat);
       
% Read the percentiles from the climatology
percentiles = ncread(climatology, 'percentiles');

% Get the index of each percentile in the climatology
for i = 1:length(selected_percentiles)
    percent_indx(i) = find(selected_percentiles(i) == percentiles);
end

% Loop over the forecast times
for lead = 0:maxlead
    % Calculate the forecasted month
    calc_mnth = tme_abs(lead + 1, 2);
    
    % Loop over all files in the current directory
    for j = 1:nfles
        % Construct the full filename
        fnme_in           = [fles(j).folder, '/', fles(j).name];
        % Load the data for the current month
        tmp               = ncread(fnme_in, varnme, [1 1 lead + 1], [Inf Inf 1]);
        % Remove all singleton dimensions
        tmp               = squeeze(tmp);
        % Save the absolute values
        X_abs(j, :, :)    = tmp;
    end
    
    % Load the mean of the current month and lead from the climatology
    climat_mn = ncread(climatology, 'mean', [1 1 calc_mnth lead + 1], [Inf Inf 1 1]);

    for j = 2:length(percent_indx) 

          % Get the threshold for the current month, lead, and percentile
          indx_up    = [1 1 calc_mnth percent_indx(j) lead + 1];
          indx_low   = [1 1 calc_mnth percent_indx(j - 1) lead + 1];
          str        = [Inf Inf 1 1 1];
          
          thresh_up  = ncread(climatology, 'thresh', indx_up, str);
          thresh_up  = squeeze(thresh_up);
          thresh_up  = reshape(thresh_up, [1, size(thresh_up, 1), ...
                                                      size(thresh_up, 2)]);
                                                     
          
          thresh_low = ncread(climatology, 'thresh', indx_low, str);
          thresh_low = squeeze(thresh_low);
          thresh_low = reshape(thresh_low, [1, size(thresh_low, 1), ...
                                                     size(thresh_low, 2)]);
          
                                                 
          mask = NaN(size(X_abs));
          
          if j == 2
              % WE DO SET VALUES <= thresh_low TO 1 AS WELL!
              mask(thresh_up >= X_abs) = 1;
          elseif j == length(percent_indx) 
              % WE DO SET VALUES > thresh_up TO 1 AS WELL!
              mask(thresh_low < X_abs) = 1;
          else
              mask(thresh_low < X_abs & thresh_up >= X_abs)  = 1;
          end
          
          
          
          categ_mn(lead + 1, j - 1, :, :)  = nanmean(X_abs.*mask, 1);
          categ_std(lead + 1, j - 1, :, :) = nanstd(X_abs.*mask, 1);
          P(lead + 1, j - 1, :, :)         = nansum(mask, 1)/nfles*100;
          
    end
      

    
    % Calculate the ensemble output statistics
    % 1. Ensemble mean
    ens_mn(lead + 1, :, :)      = mean(X_abs);
    % 2. Ensemble standard deviation
    ens_std(lead + 1, :, :)     = std(X_abs);
    % 3. Ensemble spread
    ens_spread(lead + 1, :, :)  = max(X_abs, [], 1) - min(X_abs, [], 1);
    % 4. Difference w.r.t. the climatology
    diff_climat(lead + 1, :, :) = squeeze(ens_mn(lead + 1, :, :)) - ...
                                                                 climat_mn;
    % 5. Climatology
    climat(lead + 1, :, :)      = climat_mn;

  
    
end


% Compute the event with the maximum probability
[max_probab, max_indx] = max(P, [], 2);

% Remove the singleton second dimension
max_probab = squeeze(max_probab);
max_indx   = squeeze(max_indx);
  


flgs_max(max_probab < 35)                    = 1;  
flgs_max(max_probab >= 35 & max_probab < 50) = ...
                     max_indx(max_probab >= 35 & max_probab < 50) * 10 + 1;
flgs_max(max_probab >= 50 & max_probab < 75) = ...
                     max_indx(max_probab >= 50 & max_probab < 75) * 10 + 2;
flgs_max(max_probab >= 75)                   = ...
                                       max_indx(max_probab >= 75) * 10 + 3;

flgs_max(ens_mn < 0.1) = NaN;

max_probab(max_probab == 100) = max_probab(max_probab == 100) - 0.1;
max_probab_out                = max_indx + max_probab/100;

probab_bounds = [selected_percentiles(1:end-1) selected_percentiles(2:end)];

nccreate(fnme_out, 'time', 'Dimensions', {'time', maxlead + 1});
nccreate(fnme_out, 'categories', 'Dimensions', {'categories', ncategs});
nccreate(fnme_out, 'lat', 'Dimensions', {'lat', length(lat)});
nccreate(fnme_out, 'lon', 'Dimensions', {'lon', length(lon)});
nccreate(fnme_out, 'percentile_bounds', 'Dimensions', {'categories', ncategs, 'nv', 2});

nccreate(fnme_out, 'probab', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', ncategs, ...
                                          'time', maxlead + 1}, ...
                                          'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);
                                      
nccreate(fnme_out, 'categ_mn', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', ncategs, ...
                                          'time', maxlead + 1}, ...
                                          'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);       

nccreate(fnme_out, 'categ_std', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', ncategs, ...
                                          'time', maxlead + 1}, ...
                                          'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);   
                                      
nccreate(fnme_out, 'max_probab_categ', 'Dimensions', {'lon', length(lon), ...
                                               'lat', length(lat), ...
                                               'time', maxlead + 1}, ...
                                               'Format', 'netcdf4_classic', ...
                                               'DeflateLevel', 6);
                                           
nccreate(fnme_out, 'max_probab', 'Dimensions', {'lon', length(lon), ...
                                               'lat', length(lat), ...
                                               'time', maxlead + 1}, ...
                                               'Format', 'netcdf4_classic', ...
                                               'DeflateLevel', 6);
                                           
nccreate(fnme_out, 'ens_mn', 'Dimensions', {'lon', length(lon), ...
                                            'lat', length(lat), ...
                                            'time', maxlead + 1}, ...
                                            'Format', 'netcdf4_classic', ...
                                            'DeflateLevel', 6);
                                        
                                        
nccreate(fnme_out, 'ens_std', 'Dimensions', {'lon', length(lon), ...
                                            'lat', length(lat), ...
                                            'time', maxlead + 1}, ...
                                            'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);     
                                        
nccreate(fnme_out, 'ens_spread', 'Dimensions', {'lon', length(lon), ...
                                            'lat', length(lat), ...
                                            'time', maxlead + 1}, ...
                                            'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);
                                        
                                        
nccreate(fnme_out, 'diff_climat', 'Dimensions', {'lon', length(lon), ...
                                            'lat', length(lat), ...
                                            'time', maxlead + 1}, ...
                                            'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);
                                        
nccreate(fnme_out, 'climat', 'Dimensions', {'lon', length(lon), ...
                                            'lat', length(lat), ...
                                            'time', maxlead + 1}, ...
                                            'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);
                                      

ncwriteatt(fnme_out, 'time', 'standard_name', 'time');
ncwriteatt(fnme_out, 'time', 'units', tme_unit);
ncwriteatt(fnme_out, 'lat', 'standard_name', 'latitude');
ncwriteatt(fnme_out, 'lat', 'units', 'degrees_north');
ncwriteatt(fnme_out, 'lon', 'standard_name', 'longitude');
ncwriteatt(fnme_out, 'lon', 'units', 'degrees_east');
ncwriteatt(fnme_out, 'categories', 'bounds', 'percentile_bounds');

ncwrite(fnme_out, 'lat', lat);
ncwrite(fnme_out, 'lon', lon);
ncwrite(fnme_out, 'time', tme(1:maxlead + 1));


ncwrite(fnme_out, 'categories', (1:ncategs));
ncwrite(fnme_out, 'percentile_bounds', probab_bounds);
ncwrite(fnme_out, 'probab', permute(P, [3 4 2 1]));
ncwrite(fnme_out, 'categ_mn', permute(categ_mn, [3 4 2 1]));
ncwrite(fnme_out, 'categ_std', permute(categ_std, [3 4 2 1]));
ncwrite(fnme_out, 'max_probab_categ', permute(flgs_max, [2 3 1]));
ncwrite(fnme_out, 'max_probab', permute(max_probab_out, [2 3 1]));
ncwrite(fnme_out, 'ens_mn', permute(ens_mn, [2 3 1]));
ncwrite(fnme_out, 'ens_std', permute(ens_std, [2 3 1]));
ncwrite(fnme_out, 'ens_spread', permute(ens_spread, [2 3 1]));
ncwrite(fnme_out, 'diff_climat', permute(diff_climat, [2 3 1]));
ncwrite(fnme_out, 'climat', permute(climat, [2 3 1]));



        

% Transform the absolute values in the forecasts to probabilities


