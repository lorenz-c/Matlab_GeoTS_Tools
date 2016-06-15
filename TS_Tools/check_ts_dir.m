function [adtes, mdtes, ddtes] = check_ts_dir(directory, sdte, edte, tscale, prefix, tme_frmt)

if iscell(prefix)
    files = dir([directory, prefix{1}, '*', prefix{2}]);
elseif isstr(prefix)
    files = dir([directory, prefix, '*.*']);
end

files = cellstr(char(files.name));
files(strcmp(files, '.'))  = [];
files(strcmp(files, '..')) = [];

for i = 1:length(files)
    fle_dtes{i, :} = files{i}(length(prefix{1})+1:length(prefix{1})+length(tme_frmt));
end

% Transform the dates from the files into some numbers
fle_dtes_num = datenum(fle_dtes, tme_frmt);

% Create a vector with reference dates
ref_dtes = dtevec(sdte, edte, tscale);
ref_dtes = datenum(ref_dtes);

a = 1;
m = 1;
d = 1;

adtes = {};
mdtes = {};
ddtes = {};

for i = 1:length(ref_dtes)
    indx = find(ref_dtes(i) == fle_dtes_num);
    
    if length(indx) == 1     
        adtes{a, 1} = datestr(ref_dtes(i), tme_frmt);
        a           = a + 1;
    elseif isempty(indx)    
        mdtes{m, 1} = datestr(ref_dtes(i), tme_frmt);
        m           = m + 1;
    elseif length(indx) > 1
        ddtes{d, 1} = datestr(ref_dtes(i), tme_frmt);
        d           = d + 1;
    end
end

if d == 1
    disp('No dublicate dates!')
else
    disp(['Duplicate dates: ', strjoin(ddtes)])
end

if m == 1
    disp('No missing dates!')
else
    disp(['Missing dates: ', strjoin(mdtes)])
end

        
        
    
    







