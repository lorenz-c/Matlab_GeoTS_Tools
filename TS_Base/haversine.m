function s = haversine(lat1, lon1, lat2, lon2, R, israd)
%--------------------------------------------------------------------------
% The function coputes the metric distance between two (or more) locations, 
% given by their latitude and longitude, on a sphere. 
%--------------------------------------------------------------------------
% INPUT:
% - lat1,lon1 Latitude and longitude of the first point(s)
% - lat2,lon2 Latitude and longitude of the second point(s)
% - R         Radius of the Earth (default 6378137m)
% - israd     Boolean variable which tells the function if the lats and
%             lons are in radians 
%--------------------------------------------------------------------------
% OUTPUT:
% - s         Scalar or vector with the metric distances between the two
%             (or more) locations
%--------------------------------------------------------------------------
% EXAMPLE:
% >> s = haversine(10, 20, 15, 30, 6378137, false);
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         October 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: 
%--------------------------------------------------------------------------
if nargin < 6, israd = false; end
if nargin < 5, R = 6378137; end

% Transform degrees to radians
if israd == false
    lat1 = lat1.*pi/180;
    lat2 = lat2.*pi/180;
    
    lon1 = lon1.*pi/180;
    lon2 = lon2.*pi/180;
end

delta_lat = abs(lat1 - lat2);
delta_lon = abs(lon1 - lon2);

a = sin(delta_lat/2).^2 + cos(lat1).*cos(lat2).*sin(delta_lon/2).^2;
c = 2*atan2(sqrt(a), sqrt(1-a));
s = R*c;




