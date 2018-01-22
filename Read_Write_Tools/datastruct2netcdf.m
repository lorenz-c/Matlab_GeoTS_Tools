function [] = datastruct2netcdf(inpt, outnme, compression, chunk, mval_out, ref_dte, overwrite)
% The function converts a CF-datastructure to netCDF4. It can be decided if
% all varaibles (varargin = []) or only selected variables should be
% transformed. 
% -------------------------------------------------------------------------
% Input (required):
% - inpt        CF-conform data structure
% - outnme      Name of the netCDF-file      
% - mval_out    Missing value in the netCDF-file.
% -------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
% -------------------------------------------------------------------------
% Uses: 
% -------------------------------------------------------------------------
% Updates: - 04.12.2015: Added support for time- and climatology-bounds
% -------------------------------------------------------------------------
% Set the default missing value in the output file
if nargin < 7, overwrite = true; end
if nargin < 6, ref_dte = -999; end
if nargin < 5, mval_out = 1e+20; end
if nargin < 4, chunk = []; end
if nargin < 3, compression = 8; end

warning('off', 'backtrace');

disp(['datastruct2netcdf: writing file ', outnme])
% Check if inpt contains a time-vector
if isfield(inpt.Data, 'time') 
    
    if ~isfield(inpt, 'TimeStamp')
        inpt.TimeStamp = datenum(inpt.Data.time);
    end
    
    % Number of timesteps
    nts = length(inpt.TimeStamp);
    
    % Get the first date
    if ref_dte == -999
        first_date = inpt.Data.time(1, :);
    elseif length(ref_dte) == 6
        first_date = ref_dte;
    else
        error('Unknown format of reference date')
    end
    
    if length(inpt.TimeStamp) > 1
        
        if abs(inpt.TimeStamp(2) - inpt.TimeStamp(1)) < 1
            abs(inpt.TimeStamp(2) - inpt.TimeStamp(1))
            % Hourly Data --> convert yyyy/mm/dd/hh/mm/ss into hours since 
            % first_date
            % Compute number of seconds between the date vectors
            scnds = etime(inpt.Data.time, ones(nts, 1)*first_date);
            % Transform seconds into hours
            hrs   = scnds/(60*60);
        
            % Replace the time-dimension with the new vector
            inpt.Data.time = hrs;
        
            % Replace the time-unit with the CF-conform hours since ...
            inpt.Variables.time.units = ['hours since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];

            if isfield(inpt.Variables, 'climatology_bounds')
                n_cts = size(inpt.Data.climatology_bounds, 1);
            
                bnds_left  = etime(inpt.Data.climatology_bounds(:, 1:6), ...
                                                ones(n_cts, 1)*first_date);
                bnds_right = etime(inpt.Data.climatology_bounds(:, 7:12), ...
                                                ones(n_cts, 1)*first_date);
            
                inpt.Data.climatology_bounds = ...
                                            [bnds_left bnds_right]/(60*60);
                inpt.Variables.climatology_bounds.units = ...
                                   ['hours since ', datestr(first_date, ...
                                                   'yyyy-mm-dd HH:MM:SS')];
            elseif isfield(inpt.Variables, 'time_bounds')
                n_cts = size(inpt.Data.time_bounds, 1);
            
                bnds_left  = etime(inpt.Data.time_bounds(:, 1:6), ...
                                                ones(n_cts, 1)*first_date);
                bnds_right = etime(inpt.Data.time_bounds(:, 7:12), ...
                                                ones(n_cts, 1)*first_date);
            
                inpt.Data.time_bounds = [bnds_left bnds_right]/(60*60);
                inpt.Variables.time_bounds.units = ['hours since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];
            end
            
        else
            % For daily, monthly, and annual data, convert yyyy/mm/dd/hh/mm/ss 
            % into days since first_date
            dys = daysact(ones(nts, 1)*datenum(first_date), inpt.TimeStamp);
            % Replace the time-dimension with the new vector
            inpt.Data.time = dys;
            % Replace the time-unit with the CF-conform days since ...
            inpt.Variables.time.units = ['days since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];
        
            if isfield(inpt.Variables, 'climatology_bounds')
                n_cts = size(inpt.Data.climatology_bounds, 1);
            
                bnds_left  = daysact(ones(n_cts, 1)*datenum(first_date), ...
                            datenum(inpt.Data.climatology_bounds(:, 1:6)));
                bnds_right = daysact(ones(n_cts, 1)*datenum(first_date), ...
                           datenum(inpt.Data.climatology_bounds(:, 7:12)));
           
                inpt.Data.climatology_bounds = [bnds_left bnds_right];
                inpt.Variables.climatology_bounds.units = ...
                                    ['days since ', datestr(first_date, ...
                                                   'yyyy-mm-dd HH:MM:SS')];
            elseif isfield(inpt.Variables, 'time_bounds')
                n_cts = size(inpt.Data.time_bounds, 1);
            
                bnds_left  = daysact(ones(n_cts, 1)*datenum(first_date), ...
                                   datenum(inpt.Data.time_bounds(:, 1:6)));
                bnds_right = daysact(ones(n_cts, 1)*datenum(first_date), ...
                                  datenum(inpt.Data.time_bounds(:, 7:12)));
           
                inpt.Data.time_bounds = [bnds_left bnds_right];
                inpt.Variables.time_bounds.units = ['days since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];
            end
        end
            
    else
        dys = daysact(ones(nts, 1)*datenum(first_date), inpt.TimeStamp);
        % Replace the time-dimension with the new vector
        inpt.Data.time = dys;
        % Replace the time-unit with the CF-conform days since ...
        inpt.Variables.time.units = ['days since ', datestr(first_date, ...
                                                   'yyyy-mm-dd HH:MM:SS')];
        
        if isfield(inpt.Variables, 'climatology_bounds')
            n_cts = size(inpt.Data.climatology_bounds, 1);
            
            bnds_left  = daysact(ones(n_cts, 1)*inpt.TimeStamp(1), ...
                            datenum(inpt.Data.climatology_bounds(:, 1:6)));
            bnds_right = daysact(ones(n_cts, 1)*inpt.TimeStamp(1), ...
                           datenum(inpt.Data.climatology_bounds(:, 7:12)));
            
            inpt.Data.climatology_bounds = [bnds_left bnds_right];
            inpt.Variables.climatology_bounds.units = ['days since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];
        elseif isfield(inpt.Variables, 'time_bounds')
            n_cts = size(inpt.Data.time_bounds, 1);
            
            bnds_left  = daysact(ones(n_cts, 1)*inpt.TimeStamp(1), ...
                                   datenum(inpt.Data.time_bounds(:, 1:6)));
            bnds_right = daysact(ones(n_cts, 1)*inpt.TimeStamp(1), ...
                                  datenum(inpt.Data.time_bounds(:, 7:12)));
            
            inpt.Data.time_bounds = [bnds_left bnds_right];
            inpt.Variables.time_bounds.units = ['days since ', ...
                               datestr(first_date, 'yyyy-mm-dd HH:MM:SS')];
        end
    end
   
end


%--------------------------------------------------------------------------
%                               DIMENSIONS 
%--------------------------------------------------------------------------
% Read the dimensions in the datastructure
dims = fieldnames(inpt.Dimensions);

% Number of dimensions
nr_dims = length(dims);

% Get some MetaData
MetaData = inpt.DataInfo;
          
% Create a netcdf-file
if exist(outnme, 'file') > 0
    if overwrite == true
        warning(['datastruct2netcdf: Found file ', outnme, ' in the directory. Overwriting...'])
        delete(outnme)
    elseif overwrite == false
        error(['datastruct2netcdf: Found file ', outnme, ' in the directory. Exit..'])
    end
end

mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create(outnme, mode);

% First, get the variable ID for the GLOBAL Attributes
varid_glbl = netcdf.getConstant('GLOBAL');

% Get the names of all global attributes in the input structure
glbl_Atts = fieldnames(MetaData);

% Update the file history
% new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
%                                '; MATLAB TS-Tools: datastruct2netcdf.m'];    
%        
% MetaData.history = sprintf([new_hist, ' \n', MetaData.history]);
 
% Copy the names and the corresponding values to the output file
for i = 1:length(glbl_Atts)
    netcdf.putAtt(ncid, varid_glbl, glbl_Atts{i, :}, ...
                                                  MetaData.(glbl_Atts{i}));
end

% Give each dimension in the input file a corresponding ID in the
% netCDF-file
for i = 1:nr_dims
    if isinf(inpt.Dimensions.(dims{i}))
        dim_id(i) = netcdf.defDim(ncid, dims{i}, ...
                                       netcdf.getConstant('NC_UNLIMITED'));                                 
    else
        dim_id(i) = netcdf.defDim(ncid, dims{i}, ...
                                                inpt.Dimensions.(dims{i}));  
    end
end

%--------------------------------------------------------------------------
%                               VARIABLES 
%--------------------------------------------------------------------------

% Read the variables in the datastructure
vars = fieldnames(inpt.Variables);

% Let's look at the differet variables
for i = 1:length(vars)
    % First, check if an identifier for missing values is present
    if isfield(inpt.Variables.(vars{i}), 'FillValue')
        mval = inpt.Variables.(vars{i}).FillValue;
    elseif isfield(inpt.Variables.(vars{i}), 'missingValue')
        mval = inpt.Variables.(vars{i}).missingValue;
    else
        warning(['datastruct2netcdf: No missing value attribute for variable ', vars{i}])
        mval = [];
    end
    
    % Now, copy all data from the input file to the variable bigmat; if the
    % variable does not have any data, it is assumed to be a
    % "container"-variable, which just contains some metadata.
    if isfield(inpt.Data, (vars{i}))
        if iscell(inpt.Data.(vars{i}))
            if ischar(inpt.Data.(vars{i}){1})
                bigmat{i} = char(inpt.Data.(vars{i}));
            end
        elseif isnumeric(inpt.Data.(vars{i}))
            bigmat{i} = double(inpt.Data.(vars{i}));
        end
    else
        bigmat{i} = [];
    end
    
    % Change the missing values in bigmat to 1e+20
    if ~isempty(mval)
        if ~isnan(mval)
            bigmat{i}(bigmat{i} == mval) = mval_out;
        elseif isnan(mval)
            bigmat{i}(isnan(bigmat{i}))  = mval_out;
        end
    end
    
    % Read the dimensions for each variable
    if isfield(inpt.Variables.(vars{i}), 'dimensions')
        data_dims = inpt.Variables.(vars{i}).dimensions;
    
        % Now, get the dimension IDs for the current variable
        if ~isempty(data_dims)
            for j = 1:length(data_dims)
                data_dims{j}
                var_dim_ids{i}(j) = find(ismember(dims, data_dims{j}));
            end
        
            % Get the dimension IDs in the netCDF-file
            var_dim_ids{i} = dim_id(var_dim_ids{i});
    
        
            % TWEAK: The ID which should appear first in the file should be 
            % written as the last dimension:
            var_dim_ids{i} = fliplr(var_dim_ids{i});
                                                     
        else
            var_dim_ids{i} = [];
        end
    else
        var_dim_ids{i} = [];
    end
        
    
    % Create a variable ID for each variable
    % TBA: Allow different precission for different variables!
    if isfield(inpt.Variables.(vars{i}), 'nctype')
        data_var_id(i) = netcdf.defVar(ncid, vars{i}, ...
                         inpt.Variables.(vars{i}).nctype, var_dim_ids{i});  
    else
        data_var_id(i) = netcdf.defVar(ncid, vars{i}, 'NC_DOUBLE', ...
                                                           var_dim_ids{i});
    end
                                                       
    % Set the netCDF4-compression level (values 0 - 9; 0 -> no compression, 
    % 9 -> most compression (takes longer))                                         
    netcdf.defVarDeflate(ncid, data_var_id(i), true, true, compression);
    
    if ~isempty(chunk)
       if isfield(chunk, vars{i})
            netcdf.defVarChunking(ncid, data_var_id(i), 'CHUNKED', ...
                                                          chunk.(vars{i}));
        end
    end
    
    % Get the variable attributes for the current variable
    data_Atts = fieldnames(inpt.Variables.(vars{i}));
    
    % Copy the different attributes to the netCDF-file. However, there
    % are some exceptions (e.g. FillValue, time, ...).
    for j = 1:length(data_Atts)
        if ~strcmp(data_Atts{j}, 'dimensions') & ~strcmp(data_Atts{j}, 'nctype')
            if strcmp(data_Atts{j}, 'FillValue')
                netcdf.defVarFill(ncid, data_var_id(i), false, mval_out);
     %       elseif strcmp(vars{i}, 'time') & strcmp(data_Atts{j}, 'units')
     %            netcdf.putAtt(ncid,  data_var_id(i), 'units', ...
     %          ['days since ', datestr(first_date, 'yyyy-mm-dd HH:MM:SS')]);      
            else
                inpt.Variables.(vars{i}).(data_Atts{j})
                netcdf.putAtt(ncid, data_var_id(i), data_Atts{j}, ...
                                  inpt.Variables.(vars{i}).(data_Atts{j}));
            end
        end
    end
end    

% End definitions in the netCDF-file 
netcdf.endDef(ncid);

%--------------------------------------------------------------------------
%                                   DATA
%--------------------------------------------------------------------------

% Check if the input file contains a time-dimension
tme_dim_id = find(ismember(dims, 'time')); 

if ~isempty(tme_dim_id)
    tme_var_id = find(ismember(vars, 'time'));
    netcdf.putVar(ncid, data_var_id(tme_var_id), 0, ...
                           length(bigmat{tme_var_id}), bigmat{tme_var_id});
end

for i = 1:length(vars)
    disp(['Writing variable ', vars{i}])
    if ~isempty(bigmat{i})
        if ~strcmp(vars{i}, 'time')
            if isfield(inpt.Variables.(vars{i}), 'nctype')
                if strcmp(inpt.Variables.(vars{i}).nctype, 'NC_CHAR')   
                    netcdf.putVar(ncid, data_var_id(i), char(bigmat{i}'));
                else
                    if length(var_dim_ids{i}) > 1
                        perm = 1:length(var_dim_ids{i});
                        netcdf.putVar(ncid, data_var_id(i), ...
                                         permute(bigmat{i}, fliplr(perm)));                                                 
                    else
                        netcdf.putVar(ncid, data_var_id(i), bigmat{i});
                    end
                end
            else
                if length(var_dim_ids{i}) > 1
                    perm = 1:length(var_dim_ids{i});
                    netcdf.putVar(ncid, data_var_id(i), ...
                                         permute(bigmat{i}, fliplr(perm)));   
                else
                    netcdf.putVar(ncid, data_var_id(i), bigmat{i});
                end
            end
        end
    end
end

netcdf.close(ncid);


















