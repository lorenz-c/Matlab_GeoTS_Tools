function A = area_wghts(theta, lambda, otpt, method, R, israd)
%--------------------------------------------------------------------------
% The function computes the area of grid cells according to their co-
% latitude and longitude. Therefore, four different methods can be selected
% and it can be chosen if the output should be formated as a vector or
% matrix. 
%--------------------------------------------------------------------------
% INPUT:
% - theta     Co-latitude of the pixel center (deg)
% - lambda    Longitude of the pixels (deg)
% - method    Method for computing the pixel area. Can be set to regular 
%             (default), cos, haversine, or vincenty
% - R         Radius of the Earth (default: 6378137m)
% - israd     Boolean variable which tells the function if the provided 
%             latitudes and longitudes are in degrees (default) or radians.
%--------------------------------------------------------------------------
% OUTPUT:
% - A         Area of the grid cells. The unit depends on the unit of the 
%             Earth's radius. If no radius is provided, the area is given 
%             in m^2 by default.  
%--------------------------------------------------------------------------
% EXAMPLE:
% >> A = area_wghts((89.75:-0.5:-89.75)', 179.75:0.5:179.75, 'regular'); 
%--------------------------------------------------------------------------
% Author:       Christof Lorenz (IMK-IFU)
% Date:         December 2015
% Collection:   Matlab TS-Tools 
% Version:      0.1
%--------------------------------------------------------------------------
% Uses: haversine.m, vincenty.m
%--------------------------------------------------------------------------

if nargin < 6, israd = false; end
if nargin < 5, R = 6378137; end
if nargin < 4, method = 'regular'; end
if nargin < 3, otpt = 'vec'; end

% Make sure that the thetas (lambdas) are arranged in a column (row) vector
if isrow(theta), theta = theta'; end
if iscolumn(lambda), lambda = lambda'; end

% Compute the angular sidelength of the pixels
dlat = abs(theta(2) - theta(1));
dlon = abs(lambda(2) - lambda(1));

% Compute the latitudes of the top and bottom edges of each pixel
lat1 = (theta + dlat/2);
lat2 = (theta - dlat/2);

% Compute the longitudes of the left and right edges of each pixel
lon1 = (lambda - dlon/2);
lon2 = (lambda + dlon/2);

% Transform the coordinates from degrees to radians
if israd == false
    rho = pi/180;
    
    dlat = dlat*rho;
    dlon = dlon*rho;
    
    lat1 = lat1*rho;
    lat2 = lat2*rho;
    
    lon1 = lon1*rho;
    lon2 = lon2*rho;
end

if strcmp(method, 'haversine')
    
    for i = 1:length(theta)
        sn(i,1) = haversine(lat1(i), lon1(i), lat1(i), lon2(i), R, true);
        ss(i,1) = haversine(lat2(i), lon1(i), lat2(i), lon2(i), R, true);
        sw(i,1) = haversine(lat1(i), lon1(i), lat2(i), lon1(i), R, true);
        se(i,1) = haversine(lat1(i), lon2(i), lat2(i), lon2(i), R, true);
    end
    
    A = (1/2*(sn + ss)).*(1/2*(se + sw));

elseif strcmp(method, 'regular')

    A = abs(dlon*R^2.*(sin(lat1) - sin(lat2)));
    
elseif strcmp(method, 'cos')

    A = cos(theta*pi/180)';
    
elseif strcmp(method, 'vincenty')
    
    for i = 1:length(theta)
        sn(i,1)  = vincenty(lat1(i), lat1(i), dlon*rho);
        ss(i,1)  = vincenty(lat2(i), lat2(i), dlon*rho);
        sew(i,1) = vincenty(lat1(i), lat2(i), 0);
    end

    % Interpolation for equatorial values
    snnan = find(isnan(sn));
    ssnan = find(isnan(ss));
    
    sn(snnan) = (sn(snnan-1) + sn(snnan+1))/2;
    ss(ssnan) = (ss(ssnan-1) + ss(ssnan+1))/2;
    
    A = 1/2*(sn + ss).*sew;
end

if strcmp(otpt, 'mat')
    if isrow(A)
        A = A';
    end
    A = A*ones(1, length(lambda));
end

    


  