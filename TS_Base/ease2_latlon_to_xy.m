function [X, Y] = ease2_latlon_to_xy(lat, lon)


% Parameter for the EASE-Grid 2.0

% Equatorial radius (WGS 84) 
a = 6378137.0; 

% Inverse Flattening
f = 1/298.257223563;

% Map Reference Latitude
phi0 = 0;

% Map reference Longitude
lam0 = 0;

% Eccentricity of the ellipsoid
e = sqrt(2*f - f^2);

% Standard parallel
phi1 = 30*pi/180;

phi = lat*pi/180;
lam = lon*pi/180;


k0 = cos(phi1)/sqrt(1 - e^2*sin(phi1)^2);

q_phi = (1 - e^2)*((sin(phi)./(1 - e^2*sin(phi).^2)) - 1/(2*e)* ...
    log((1 - e*sin(phi))./(1 + e*sin(phi))));



X = a*k0*(lam - lam0);
Y = a*q_phi./(2*k0);



