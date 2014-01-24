function [RE] = Z_RE(lat1,lat2)
% 11/28/2012 Parker MacCready
%
% [RE] = Z_RE(lat1,lat2)
%
% calculate the Earth radius (m) at the average latitude of the domain
% (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid
%
% INPUT: lat1 and lat2 (degrees) [or just lat1]
%
% OUTPUT: Earth radius (m) at the mean latitute

if nargin==1
    latm = lat1;
elseif nargin==2
    latm=(lat1+lat2)/2;
end

a = 6378.137 * 1000; % equatorial radius (m)
b = 6356.7523 * 1000; % polar radius (m)
cl=cos(pi*latm/180);
sl=sin(pi*latm/180);
RE = sqrt(((a*a*cl)^2 + (b*b*sl)^2) / ((a*cl)^2 + (b*sl)^2));
