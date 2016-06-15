function fld_out = structcat(dim, varargin)

if iscell(varargin{1,1})
    for i = 1:length(varargin{1,1})
        tmp{i} = varargin{1}{i};
    end
    varargin = tmp;
end


fld_out = varargin{1};

if strcmp(dim, 'time')
    for i = 1:length(varargin)
        istime            = istimevar(varargin{i});
        vars              = fieldnames(varargin{i}.Variables);
        vars(istime ~= 1) = [];
        
        tme_pos = getdimpos(varargin{i}, vars, 'time');
       
        for j = 1:length(vars)
            if i == 1
                fld_out.Data.(vars{j}) = varargin{1}.Data.(vars{j});
            else
                if tme_pos(j) == 1
                    fld_out.Data.(vars{j}) = [fld_out.Data.(vars{j}); ...
                                              varargin{i}.Data.(vars{j})];
                elseif tme_pos(j) == 2
                    fld_out.Data.(vars{j}) = [fld_out.Data.(vars{j}) ...
                                              varargin{i}.Data.(vars{j})];
                end
            end
        end
    end
    fld_out.TimeStamp = datenum(fld_out.Data.time);    
    if isinf(varargin{1}.Dimensions.time)
        fld_out.Dimensions.time = Inf;
    else
        fld_out.Dimensions.time = length(fld_out.TimeStamp);
    end
end
        
    

