function ddt = struct_cntrl_diff(inpt, vars, dt)
% The function computes the first derivative w.r.t. time of a given input
% dataset through the method of central differences (see, e.g., 
% http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture4.pdf).
% Therefore, it is checked if the inpt-data has continuous time-steps. The
% derivative of the first and last time-step is computed through forward
% and backward differences, respectively.
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Datastructure from which the derivative should be computed
%--------------------------------------------------------------------------
% OUTPUT:
% - ddt         Matlab datastructure which contains the first derivative
%               for each variable in the vars-list
% - vars        Cell-array with variables, from which the derivative should
%               be computed. If set to 'all' (or empty), the function will 
%               compute the derivative for all non-fixed variables.
% - dt          Lenght of the time-period between two dates. By default,
%               dt is set to 1, i.e. the unit of the output corresponds
%               directly to the temp. resolution of the inpt (per hour, 
%               day, month, ...)
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         January
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: isfixedvar.m
%--------------------------------------------------------------------------

if nargin < 2, vars = 'all'; end
if nargin < 3, dt = 1; end

% Create a new datastructure with the metadata and dimensions from the
% input data
ddt.DataInfo   = inpt.DataInfo;
ddt.Dimensions = inpt.Dimensions;

if strcmp(vars, 'all')
    vars = fieldnames(inpt.Variables);
    
    for i = 1:length(vars)
        isfixed(i) = isfixedvar(vars{i});
    end
    % Save the fixed variables
    vars_fixed = vars(isfixed == 1);
    % Remove the fixed variables from the variables-list for the following
    % computations
    vars(isfixed == 1) = [];
end

% Check the length of the dataset
inpt_tmp = trunc_TS(inpt, inpt.Data.time(1, :), inpt.Data.time(end, :));

% Copy the fixed variables to the new datastructure
ddt = copyvars(ddt, inpt_tmp, vars_fixed);


                                               
for i = 1:length(vars)
    % Get the "position" of the time-dimension
    time_id(i) = find(ismember(inpt.Variables.(vars{i}).dimensions, ...
                                                               'time'), 1);
    % Get the number of dimensions for each variable
    nr_dims = length(inpt.Variables.(vars{i}).dimensions);
    
    if nr_dims <= 2
        if time_id(i) == 2
            fld = inpt_tmp.Data.(vars{i})';
        else
            fld = inpt_tmp.Data.(vars{i});
        end
        
        % Add the first and last time-step to the beginning and end of the
        % time-series
        fld_enh = [fld(1, :); fld; fld(end, :)];
        
        % Compute the derivative
        ddt.Data.(vars{i}) = 1/(2*dt)*(fld_enh(3:end, :) - ...
                                                      fld_enh(1:end-2, :));
        
        
        
    elseif nr_dims == 3
        if time_id(i) == 1
            fld = inpt_tmp.Data.(vars{i});
        else
            error('First dimension must be "time"!')
        end
        
        % Add the first and last time-step to the beginning and end of the
        % time-series
        fld_enh = [fld(1, :, :); fld; fld(end, :, :)];
        
        % Compute the derivative
        ddt.Data.(vars{i}) = 1/(2*dt)*(fld_enh(3:end, :, :) - ...
                                                   fld_enh(1:end-2, :, :));
        
    elseif nr_dims == 4
        if time_id(i) == 1
            fld = inpt_tmp.Data.(vars{i});
        else
            error('First dimension must be "time"!')
        end
        
        % Add the first and last time-step to the beginning and end of the
        % time-series
        fld_enh = [fld(1, :, :, :); fld; fld(end, :, :, :)];
        
        % Compute the derivative
        ddt.Data.(vars{i}) = 1/(2*dt)*(fld_enh(3:end, :, :, :) - ...
                                                fld_enh(1:end-2, :, :, :));
    end
    
    % Copy the variable's metadata to the output
    ddt.Variables.(vars{i}) = inpt.Variables.(vars{i});
    
    % Update the source-attribute of the variables
    new_srce = 'First derivative w.r.t. time using central differences.';
    
    if isfield(ddt.Variables.(vars{i}), 'source')
        ddt.Variables.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                          ddt.Variables.(vars{i}).source]);
    else
        ddt.Variables.(vars{i}).source = new_srce;
    end
end

% Copy the time-data to the output
ddt = copyvars(ddt, inpt_tmp, {'time'});

% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                 '; MATLAB TS-Tools: struct_cntrl_diff.m'];

if isfield(ddt.DataInfo, 'history')           
    ddt.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                    ddt.DataInfo.history]);
else
    ddt.DataInfo.history = new_hist;
end        

warning('off', 'backtrace')
warning('Units have changed depending on the temporal step size!')
warning('Please check manually!')
    