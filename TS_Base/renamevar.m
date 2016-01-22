function out = renamevar(inpt, var_old, var_new)

% Copy the inpt to the otpt
out = inpt;

if isstr(var_old)
    var_old = {var_old};
    var_new = {var_new};
end

vars = fieldnames(out.Variables);

for i = 1:length(var_old)
    % Check if var_old is used as a dimension
    if isfield(out.Dimensions, var_old{i})
        for j = 1:length(vars)
            var_dims = out.Variables.(vars{j}).dimensions;
            var_indx = find(ismember(var_dims, var_old{i}) == 1);
        
            if ~isempty(var_indx)
                var_dims{var_indx} = var_new{i};
                out.Variables.(vars{j}).dimensions = var_dims;
            end
        end   

        % Copy the size of the dimension to the new variable name
        out.Dimensions.(var_new{i}) = out.Dimensions.(var_old{i});
        % Remove the old variable
        out.Dimensions           = rmfield(out.Dimensions, var_old{i});
    end

    out.Variables.(var_new{i}) = out.Variables.(var_old{i});
    out.Variables           = rmfield(out.Variables, var_old{i});

    out.Data.(var_new{i})      = out.Data.(var_old{i});
    out.Data                = rmfield(out.Data, var_old{i});
    
    % Update the variable list
    vars = fieldnames(out.Variables);
end


new_hist = [datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                                         '; MATLAB TS-Tools: renamevar.m'];    
        
out.DataInfo.history = sprintf([new_hist, ' \n', out.DataInfo.history]);