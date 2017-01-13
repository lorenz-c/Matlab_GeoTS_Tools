function out = struct_movav(inpt, vars, window)


if nargin < 2, vars = 'all'; end
if nargin < 3, window = 3; end

% Create a new datastructure with the metadata and dimensions from the
% input data
out.DataInfo   = inpt.DataInfo;
out.Dimensions = inpt.Dimensions;

vars = fieldnames(inpt.Variables);
    
% Get all fixed variables
for i = 1:length(vars)
    isfixed(i) = isfixedvar(vars{i});
end

% Save the fixed variables
out = copyvars(out, inpt, vars(isfixed == 1));

% Remove the fixed variables from the variables-list for the following
% computations
vars(isfixed == 1) = [];
    
% Get all variables which have a 'time'-dimension
istime = istimevar(inpt, vars);

% Save the variables without the time-dimension
out = copyvars(out, inpt, vars(istime == 0));

% Remove the variables without the time-dimension
vars(istime == 0) = [];

% Get the "position" of the time-dimension  
timepos = getdimpos(inpt, vars, 'time');

% Get the number of dimensions for each remaining variable
nrdims = getnrdims(inpt, vars);


n  = window;

for i = 1:length(vars)

    if nrdims(i) <= 2
        if timepos(i) == 2
            fld = inpt_tmp.Data.(vars{i})';
        else
            fld = inpt_tmp.Data.(vars{i});
        end
        
        fld_out = NaN(size(fld));
        
        % Depending on the length of the window, enlarge the data with the
        % first and last "row" so that the size of the output matches the
        % size of the input.
        frst = repmat(fld(1, :),   floor(n/2), 1);
        last = repmat(fld(end, :), floor(n/2), 1);

        fld_enh = [frst; fld; last];

        for j = floor(n/2) + 1:size(fld_enh, 1) - floor(n/2)
            fld_out(j-floor(n/2), :, :) = ...
                         nanmean(fld_enh(j-floor(n/2):j+floor(n/2), :), 1);
        end
        
        if timepos(i) == 2
            fld_out = fld_out';
        end
        
        
    elseif nrdims == 3
        
        if timepos(i) == 1
            fld = inpt.Data.(vars{i});
        else
            error('First dimension must be "time"!')
        end
        
        fld_out = NaN(size(fld));
        
        % Depending on the window size, add the first and last time-steps 
        % to the beginning and end of the time-series
        frst = repmat(fld(1, :, :), [floor(n/2), 1, 1]);
        last = repmat(fld(end, :, :), [floor(n/2), 1, 1]);
        
        fld_enh = [frst; fld; last];

        for j = floor(n/2) + 1:size(fld_enh, 1) - floor(n/2)
            fld_out(j-floor(n/2), :, :) = ...
                      nanmean(fld_enh(j-floor(n/2):j+floor(n/2), :, :), 1);
        end
        
    elseif nrdims == 4
        
        if timepos(i) == 1
            fld = inpt_tmp.Data.(vars{i});
        else
            error('First dimension must be "time"!')
        end
        
        fld_out = NaN(size(fld));
        
        % Depending on the window size, add the first and last time-steps 
        % to the beginning and end of the time-series
        frst = repmat(fld(1, :, :, :), [floor(n/2), 1, 1, 1]);
        last = repmat(fld(end, :, :, :), [floor(n/2), 1, 1, 1]);
        
        fld_enh = [frst; fld; last];

        for j = floor(n/2) + 1:size(fld_enh, 1) - floor(n/2)
            fld_out(j-floor(n/2), :, :, :) = ...
                   nanmean(fld_enh(j-floor(n/2):j+floor(n/2), :, :, :), 1);
        end
    end
   
    % Copy the variable's metadata to the output
    out.Variables.(vars{i}) = inpt.Variables.(vars{i});
    
    % Update the source-attribute of the variables
    new_srce = ['Moving average over time with n = ', num2str(n)];
    
    if isfield(out.Variables.(vars{i}), 'source')
        out.Variables.(vars{i}).source = sprintf([new_srce, ' \n', ...
                                          out.Variables.(vars{i}).source]);
    else
        out.Variables.(vars{i}).source = new_srce;
    end
    
    out.Data.(vars{i}) = fld_out;
    clear fld_out
end



% Update the history
new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                 '; MATLAB TS-Tools: struct_movav.m'];
if isfield(out.DataInfo, 'history')           
    out.DataInfo.history = sprintf([new_hist, ' \n', ...
                                                   out.DataInfo.history]);
else
    out.DataInfo.history = new_hist;
end        

    