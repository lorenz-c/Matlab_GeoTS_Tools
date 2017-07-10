function [lon, lat] = ease_inv(x, y)





% Parameter of WGS84
a = 6378137;
f = 1/298.257223563;

b = a * (1 - f);
e = sqrt((a^2 - b^2)/a^2);

% Parameter of the EASE 2.0 Projection

lam0 = 0;
phi1 = 30*pi/180;
phi0 = 90*pi/180;


k0 = cos(phi1)/sqrt(1-e^2*sin(phi1)^2);
qp = (1-e^2) * (sin(phi0)/(1-e^2*sin(phi0)^2) - 1/(2*e) * ...
    log((1 - e*sin(phi0))/(1 + e*sin(phi0))));

beta = asin(2*y*k0./(a*qp));

lon = lam0 + x/(a*k0);
lat = beta + (e^2/3 + 31*e^4/180 + 517*e^6/5040)*sin(2*beta) + ...
    (23*e^4/360 + 251*e^6/3780)*sin(4*beta) + ...
    (761*e^6/45360)*sin(6*beta);

lon = lon*180/pi;
lat = lat*180/pi;

