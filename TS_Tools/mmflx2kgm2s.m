function out = mmflx2kgm2s(inpt)

% The function converts a datastructure with variables in units of
% [mm/month] or [mm/day] to (COARDS and CF-conform) [kg/m2/s]
vars = fieldnames(inpt.Variables)

dtes  = datevec(inpt.TimeStamp);
yrs   = dtes(:, 1);
mnths = dtes(:, 2);

out = inpt;

for i = 1:length(vars)
    var_unit = inpt.Variables.(vars{i}).units;
    
    if strcmp(var_unit, 'mm/month')
        nrd = eomday(yrs, mnths);
        
        sze_dta    = size(inpt.Data.(vars{i}));
        sze_dta(1) = 1;
        
        nrd = repmat(nrd, sze_dta);
        
        Data_out = inpt.Data.(vars{i})./(nrd*24*60*60);
        out.Variables.(vars{i}).units = 'kg/m^2/s';
                out.Data.(vars{i}) = Data_out;
                
    elseif strcmp(var_unit, 'mm/day')
        
        Data_out = inpt.Data.(vars{i})/(24*60*60);
        out.Variables.(vars{i}).units = 'kg/m^2/s';
        out.Data.(vars{i}) = Data_out;
        
    elseif strcmp(var_unit, 'mm/hr') | strcmp(var_unit, 'mm/hour')
        
        Data_out = inpt.Data.(vars{i})/(60*60);
        out.Variables.(vars{i}).units = 'kg/m^2/s';
        out.Data.(vars{i}) = Data_out;
        
    end
    
    
    
end

if isfield(out.DataInfo, 'history')
    out.DataInfo.history = sprintf([[datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                       '; mmflx2kgm2s.m: Unit conversion'], ' \n', ...
                                        out.DataInfo.history]);
elseif isfield(out.DataInfo, 'History')
    out.DataInfo.History = sprintf([[datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                       '; mmflx2kgm2s.m: Unit conversion'], ' \n', ...
                                        out.DataInfo.History]);
else
    out.DataInfo.history = sprintf([datestr(now, 'ddd mmm dd HH:MM:SS yyyy'), ...
                       '; mmflx2kgm2s.m: Unit conversion']);
end
                
 
            
        





    
