function P = probabilistic_forecast(inpt_dir, varnme, distr, thresh, maxlead, fnme_out, lead_climat)
   
if nargin < 7
    lead_climat = 0;
end

if nargin < 5
    maxlead = 6;
end

if nargin < 4
    thresh = [0 0.1 0.33 0.66 0.9 1]';
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
ncategs = length(thresh) - 1;

% Create empty matrices for all output variables
probab      = NaN(maxlead + 1, ncategs, nlat, nlon);
categ_mn    = NaN(maxlead + 1, ncategs, nlat, nlon);
categ_std   = NaN(maxlead + 1, ncategs, nlat, nlon);
ens_mn      = NaN(maxlead + 1,          nlat, nlon);
ens_std     = NaN(maxlead + 1,          nlat, nlon);
ens_spread  = NaN(maxlead + 1,          nlat, nlon);
diff_climat = NaN(maxlead + 1,          nlat, nlon);
climat      = NaN(maxlead + 1,          nlat, nlon);
        
        
% Loop over the forecast times
for i = 1:maxlead + 1
    % Identify the actual month
    calc_mnth = tme_abs(i, 2);
    % Get the alpha-and beta parameter and climatology for the current 
    % month 
    if lead_climat == 0
        alpha     = squeeze(distr.Data.alpha(calc_mnth, :, :));
        beta      = squeeze(distr.Data.beta(calc_mnth, :, :));
        mn_climat = squeeze(distr.Data.mean(calc_mnth, :, :));
    elseif lead_climat == 1
        alpha     = squeeze(distr.Data.alpha(i, calc_mnth, :, :));
        beta      = squeeze(distr.Data.beta(i, calc_mnth, :, :));
        mn_climat = squeeze(distr.Data.mean(i, calc_mnth, :, :));
    end
        
    
    % Loop over all files in the current directory
    for j = 1:length(fles)
        % Construct the full filename
        fnme_in           = [fles(j).folder, '/', fles(j).name];
        % Load the data for the current month
        tmp               = ncread(fnme_in, varnme, [1 1 i], [Inf Inf 1]);
        % Remove all singleton dimensions
        tmp               = squeeze(tmp)';
        % Save the absolute values
        X_abs(j, :, :)    = tmp;
        % ...and perform the probability transform using the Gamma
        % distribution and the parameter from the climatology
        X_probab(j, :, :) = cdf('gamma', tmp, alpha, beta);   
    end
    
    % Loop over the percentile thresholds
    for j = 1:length(thresh)-1
        % Create a mask for saving all values that fall within a category
        mask_thresh = NaN(size(X_probab));
        
        if j == length(thresh)-1
            % For the last cateogry, also take values that are equal to the
            % upper bound
            mask_thresh(X_probab >= thresh(j) & X_probab <= thresh(j+1)) = 1;
        else
            % For the other cases, look for values that are equal or
            % greater than the lower bound and lower than the upper bound
            mask_thresh(X_probab >= thresh(j) & X_probab < thresh(j+1)) = 1;
        end
        
        % Compute the categorical mean (i.e. the mean from all ensemble
        % member that fall within a category)
        categ_mn(i, j, :, :)  = squeeze(nanmean(X_abs.*mask_thresh, 1));
        % ...and standard deviation
        categ_std(i, j, :, :) = squeeze(nanstd(X_abs.*mask_thresh, 1));
        % Finally, compute the probability of the different categories
        P(i, j, :, :)         = squeeze(nansum(mask_thresh, 1))/length(fles)*100;
        
    end
    
        
    
    % Calculate the ensemble output statistics
    % 1. Ensemble mean
    ens_mn(i, :, :)      = mean(X_abs);
    % 2. Ensemble standard deviation
    ens_std(i, :, :)     = std(X_abs);
    % 3. Ensemble spread
    ens_spread(i, :, :)  = max(X_abs, [], 1) - min(X_abs, [], 1);
    % 4. Difference w.r.t. the climatology
    diff_climat(i, :, :) = squeeze(ens_mn(i, :, :)) - mn_climat;
    % 5. Climatology
    climat(i, :, :)      = mn_climat;


    
end

probab_bounds = [thresh(1:end-1) thresh(2:end)];

nccreate(fnme_out, 'time', 'Dimensions', {'time', maxlead + 1});
nccreate(fnme_out, 'categories', 'Dimensions', {'categories', length(thresh)-1});
nccreate(fnme_out, 'lat', 'Dimensions', {'lat', length(lat)});
nccreate(fnme_out, 'lon', 'Dimensions', {'lon', length(lon)});
nccreate(fnme_out, 'percentile_bounds', 'Dimensions', {'categories', length(thresh)-1, 'nv', 2});

nccreate(fnme_out, 'probab', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', length(thresh)-1, ...
                                          'time', maxlead + 1}, ...
                                          'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);
                                      
nccreate(fnme_out, 'categ_mn', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', length(thresh)-1, ...
                                          'time', maxlead + 1}, ...
                                          'Format', 'netcdf4_classic', ...
                                          'DeflateLevel', 6);       

nccreate(fnme_out, 'categ_std', 'Dimensions', {'lon', length(lon), ...
                                          'lat', length(lat), ...
                                          'categories', length(thresh)-1, ...
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

ncwrite(fnme_out, 'categories', (1:length(thresh)-1));
ncwrite(fnme_out, 'percentile_bounds', probab_bounds);
ncwrite(fnme_out, 'probab', permute(P, [4 3 2 1]));
ncwrite(fnme_out, 'categ_mn', permute(categ_mn, [4 3 2 1]));
ncwrite(fnme_out, 'categ_std', permute(categ_std, [4 3 2 1]));
ncwrite(fnme_out, 'ens_mn', permute(ens_mn, [3 2 1]));
ncwrite(fnme_out, 'ens_std', permute(ens_std, [3 2 1]));
ncwrite(fnme_out, 'ens_spread', permute(ens_spread, [3 2 1]));
ncwrite(fnme_out, 'diff_climat', permute(diff_climat, [3 2 1]));
ncwrite(fnme_out, 'climat', permute(climat, [3 2 1]));



        

% Transform the absolute values in the forecasts to probabilities


