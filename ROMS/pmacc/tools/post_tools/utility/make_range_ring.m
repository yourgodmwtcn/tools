function circle = make_range_ring(long,lat,radius)
%------------------------------------------------
%
% modified from m_range_ring.m file in m-map toolbox to output circle
% 
% M_RANGE_RING Creates range rings on a map
%    M_RANGE_RING(LONG,LAT,RADIUS) creates a range ring of Radius (in km)
%    centered at the position specified by LONG and LAT. 

% Rich Pawlowicz (rich@ocgy.ubc.ca) 18/Dec/1998
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.

% 6/Nov/00 - eliminate returned stuff if ';' neglected (thx to D Byrne)
%-----------------------------------------------

pi180 = pi/180;
RE = Z_RE(lat,lat);
n = 100;

c = radius*1000/RE;

rlat = lat*pi180;
rlong = long*pi180;

x = sin([0:n-1]'/(n-1)*2*pi)*c;
y = cos([0:n-1]'/(n-1)*2*pi)*c;
on = ones(n,1);

Y = (asin(on*cos(c)*sin(rlat) + (on*cos(rlat)*(sin(c)./c)).*y))/pi180;
  
X = (rlong+atan2(x.*(on*sin(c)),on*(cos(rlat)*cos(c).*c) - (on*sin(rlat)*sin(c)).*y ) )/pi180;
 
circle = [X Y]; %output circle polygon

%------------------ calculate earth's radius -----------------------
function [RE] = Z_RE(lat1,lat2)
% 10/16/2008  Parker MacCready

% calculate the Earth radius at the average latitude of the domain
% (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid
a = 6378.137 * 1000; % equatorial radius (m)
b = 6356.7523 * 1000; % polar radius (m)
latm=(lat2+lat1)/2; cl=cos(pi*latm/180); sl=sin(pi*latm/180);
RE = sqrt(((a*a*cl)^2 + (b*b*sl)^2) / ((a*cl)^2 + (b*sl)^2));
