function otpt = kgm2s2mmflx(inpt, out_unit)
% The function converts a datastructure with flux variables in units of
% [kg/m2/s] to [mm/month] or [mm/day]. Please note that for the conversion
% to [mm/month], the function requires the respective month of each
% time-step. Therefore, it is assumed that the input structure contains a
% variable "time", where the dates are stored in 6 columns 
% (YYYY MM DD HH MM SS), from which the first two contain the years and
% months.
%--------------------------------------------------------------------------
% Input (required):
% - inpt        CF-conform data structure
% - out_unit    Defines the unit of the output. Possible selections are
%               'mm/month' and 'mm/day' (default)

%
% Output:
% - out     
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% If out_unit is empty, set it to [mm/day]
if nargin < 2, out_unit = 'mm/day'; end

% Get the input parameters (for updating the history)
params_in{1} = [inputname(1), ','];
params_in{2} = out_unit;

% Copy the input to the output variable for keeping all MetaData, etc.
otpt = inpt;

% Get the variables of the input 
vars = fieldnames(inpt.Variables);


if strcmp(out_unit, 'mm/month') 
    for i = 1:length(vars)
        % Does the current variable have an "units"-attribute?
        if isfield(inpt.Variables.(vars{i}), 'units')
            % Is the current variable in units of kg/m^2/s?
            if strcmp(inpt.Variables.(vars{i}).units, 'kg/m^2/s') | ...
                  strcmp(inpt.Variables.(vars{i}).units, 'kgm^-2s^-1') | ...
                       strcmp(inpt.Variables.(vars{i}).units, 'kg m-2 s-1')
                % Short command line text message
                disp(['kgm2s2mmflx.m: Found matching variable: ', vars{i}]) 
         
                % Get the years and months from the input data
                DateTime = inpt.Data.time;
                Yrs      = DateTime(:, 1);
                Mnths    = DateTime(:, 2);
            
                % Compute the number of days for each month
                nrd      = eomday(Yrs, Mnths)*24*3600;
                
                % Get the "position" of the time-dimension
                dimpos = getdimpos(inpt, vars{i}, 'time');
            
                % Multiply each data slice with the corresponding number of
                % days (x 60 seconds x 60 minutes x 24 hours)
                if dimpos == 1
                    otpt.Data.(vars{i}) = bsxfun(@times, ...
                                                inpt.Data.(vars{i}), nrd);
                elseif dimpos == 2
                    otpt.Data.(vars{i}) = bsxfun(@times, ...
                                                inpt.Data.(vars{i}), nrd');
                end
            
                % Update the variable's MetaData
                otpt.Variables.(vars{i}).units = 'mm/month';
            end
        end
    end
    
elseif strcmp(out_unit, 'mm/day')
    
    % Conversion to [mm/day]
    for i = 1:length(vars)
        % Does the current variable have an "units"-attribute?
        if isfield(inpt.Variables.(vars{i}), 'units')
            % Is the current variable in units of kg/m^2/s?
            if strcmp(inpt.Variables.(vars{i}).units, 'kg/m^2/s') | ...
                        strcmp(inpt.Variables.(vars{i}).units, 'kgm^-s^-1')
                    
                % Short command line text message    
                disp(['kgm2s2mmflx.m: Found matching variable: ', vars{i}]) 
           
                % Multiply each data slice with 60 seconds x 60 minutes x 
                % 24 hours
                otpt.Data.(vars{i}) = inpt.Data.(vars{i})*3600*24;
            
                % Update the variable's MetaData
                otpt.Variables.(vars{i}).units = 'mm/day';
            end
        end
    end
else
    error('kgm2s2mmflx.m: Unknown output unit!')
end

% Update the file history  
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
              '; MATLAB TS-Tools: kgm2s2mmflx.m; output unit: ', out_unit];     
otpt.DataInfo.history = sprintf([new_hist, ' \n', inpt.DataInfo.history]);

         
