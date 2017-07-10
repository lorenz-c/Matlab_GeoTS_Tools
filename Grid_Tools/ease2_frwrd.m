function [x, y] = ease_frwrd(lon, lat, israd)

if nargin < 3, israd = false; end

% Parameter of WGS84
a = 6378137;
f = 1/298.257223563;

b = a * (1 - f);
e = sqrt((a^2 - b^2)/a^2);

% Parameter of the EASE 2.0 Projection

lam0 = 0;
phi1 = 30*pi/180;


if israd == false
    lon = lon*pi/180;
    lat = lat*pi/180;
end


k0 = cos(phi1)/sqrt(1-e^2*sin(phi1)^2);
q  = (1-e^2) * (sin(lat)./(1-e^2*sin(lat).^2) - 1/(2*e) * ...
    log((1 - e*sin(lat))./(1 + e*sin(lat))));

x = a*k0 *(lon - lam0);
y = a*q./(2*k0);

