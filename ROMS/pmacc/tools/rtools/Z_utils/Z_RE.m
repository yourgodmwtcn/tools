function [RE] = Z_RE(lat1,lat2)
% 10/20/2011  Parker MacCready
%
% calculate the Earth radius (m) at the average latitude of the domain
% (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid
a = 6378.137 * 1000; % equatorial radius (m)
b = 6356.7523 * 1000; % polar radius (m)
latm=(lat2+lat1)/2; cl=cos(pi*latm/180); sl=sin(pi*latm/180);
RE = sqrt(((a*a*cl)^2 + (b*b*sl)^2) / ((a*cl)^2 + (b*sl)^2));
