function otpt = copyvars(inpt, source, vars)
% The function copies variables (i.e. their metadata and data) to the 
% output. 
%--------------------------------------------------------------------------
% INPUT:
% - inpt        Datastructure, to which the variables should be added
% - source      Source datastructure, from which the variables are copied
% - vars        Cell array where each cell contains a variable name. All
%               the variables in vars are copied from source to inpt
%--------------------------------------------------------------------------
% OUTPUT:
% - otpt        Matlab datastructure with the new variables
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         November 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% Copy all fields from inpt to otpt
otpt = inpt;

% Get the dimensions of the inpt datastructure
if isfield(inpt, 'Dimensions')
    vars_inpt = fieldnames(inpt.Dimensions);
else
    vars_inpt = [];
end
    
for i = 1:length(vars)
    % Get the dimensions of the source datastructure
    if isfield(source.Variables.(vars{i}), 'dimensions')
        vardims = source.Variables.(vars{i}).dimensions;
    
        if ~isempty(vars_inpt)
            % Check if all dimensions of the source variable are also present in 
            % the inpt datastructure
            for j = 1:length(vardims)
                if ~ismember(vardims{j}, vars_inpt)
                    otpt.Dimensions.(vardims{j}) = ...
                                            source.Dimensions.(vardims{j});
                end
            end
        else
            otpt.Dimensions.(vardims{i}) = source.Dimensions.(vardims{i}); 
        end
    end
    % Copy the variable's metadata
    otpt.Variables.(vars{i}) = source.Variables.(vars{i});
    % Copy the variable's data
    if isfield(source.Data, vars{i})
        otpt.Data.(vars{i})      = source.Data.(vars{i});
    end
    
    % If the time-variable is added to the otpt, we also add the respective
    % TimeStamp
    if strcmp(vars, 'time')
        if isfield(source, 'TimeStamp')
            otpt.TimeStamp = source.TimeStamp;
        else
            otpt.TimeStamp = datenum(source.Data.time);
        end
    end
end


