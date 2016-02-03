function out = findregions(inpt, region_ids, region_var, miss_regions)

if nargin < 4, miss_regions = true; end
if nargin < 3, region_var   = 'regions'; end

% Get the region ids of the input data
regions_dta = inpt.Data.(region_var);

% Go through all region_ids and search for the indices in regions_dta
for i = 1:length(region_ids)
    tmp = find(region_ids(i) == regions_dta);
    if ~isempty(tmp)
        indx(i) = tmp;
    else
        indx(i) = NaN;
    end
end

% If we want to keep the regions which are not found in the dataset we need
% the indices of the missing regions and the regions which are included in
% inpt.
if miss_regions == true
    nan_regions = find(isnan(indx));
    val_regions = find(~isnan(indx));
end
    
% Remove the missing regions from indx
indx(isnan(indx)) = [];

% Get the variables in inpt
vars = fieldnames(inpt.Variables);

for i = 1:length(vars)
    % Check if the current variable is connected to the region_var
    isregion = ismember(region_var, inpt.Variables.(vars{i}).dimensions);
    
    if isregion == 1
        % Get the position of the region dimension
        dimpos = getdimpos(inpt, vars{i}, region_var);
        % Get the number of dimensions
        dims   = length(inpt.Variables.(vars{i}).dimensions);
        
        
        % First, check if the current variable is a cell array (e.g. region
        % names, ...)
        if iscell(inpt.Data.(vars{i}))
            if dims <= 2
                if dimpos == 1
                    if miss_regions == true
                        out.Data.(vars{i})(val_regions, 1) = ...
                                                 inpt.Data.(vars{i})(indx);
                        out.Data.(vars{i})(nan_regions, 1) = {'undefined'};
                    else
                        out.Data.(vars{i}) = inpt.Data.(vars{i})(indx);
                    end
                elseif dimpos == 2
                    if miss_regions == true
                        out.Data.(vars{i})(1, val_regions) = ...
                                                 inpt.Data.(vars{i})(indx);
                        out.Data.(vars{i})(1, nan_regions) = {'undefined'};
                    else
                        out.Data.(vars{i}) = inpt.Data.(vars{i})(indx);
                    end
                end
            else
                error([vars{i}, ' has unknown data tye.'])
            end
        % If not, we can use matrix notation for the indexing    
        else
            if dims <= 2
                if dimpos == 1
                    if miss_regions == true
                        out.Data.(vars{i})(val_regions, :) = ...
                                              inpt.Data.(vars{i})(indx, :);

                        out.Data.(vars{i})(nan_regions, :) = NaN;
                    else                                
                    	out.Data.(vars{i}) = inpt.Data.(vars{i})(indx, :);
                    end
                
                elseif dimpos == 2
                    if miss_regions == true
                        out.Data.(vars{i})(:, val_regions) = ...
                                              inpt.Data.(vars{i})(:, indx);
                        out.Data.(vars{i})(:, nan_regions) = NaN;
                    else                                
                        out.Data.(vars{i}) = inpt.Data.(vars{i})(:, indx);
                    end
                end    
            end
        end
        
        % Copy the variable's metadata
        out.Variables.(vars{i}) = inpt.Variables.(vars{i});
    end
end

% If we want to keep the missing regions, we use the region_ids as new
% region_var
if miss_regions == true
    out.Data.(region_var) = region_ids;
end

% If inpt has a time variable, we copy that to the output
if ismember('time', vars)
    out = copyvars(out, inpt, {'time'});
end

% Copy the input's metadata
out.DataInfo = inpt.DataInfo;

% Finally, adjust the length of the region dimension
out.Dimensions.(region_var) = length(indx);
 
    
