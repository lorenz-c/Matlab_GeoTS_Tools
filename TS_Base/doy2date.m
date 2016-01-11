function date = doy2date(year, doy)
%--------------------------------------------------------------------------
% Transform day of year into datetime
%--------------------------------------------------------------------------
% INPUT:
% - year      Reference year
% - doy       Day of the year
%--------------------------------------------------------------------------
% OUTPUT:
% - date      Date-time matrix with YYYY MM DD HH MM SS
%--------------------------------------------------------------------------
% EXAMPLE:
% >> date = doy2date(2000, 150);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------

% 0. Check if the year is a leap-year
n_feb  = eomday(year, 2);

% 1. Create a vector with days of the year
if n_feb == 29
    doy_ref = 1:366;
else
    doy_ref = 1:365;
end

% 2. Create a vector with daily dates
dte = dtevec([year 01], [year 12], 'daily');

% 3. Now, find the index of the respective doy
doy_indx = find(doy_ref == doy);

% 4. Use the respective (doy_indx) element from the date vector
date = dte(doy_indx, :);
         
        

