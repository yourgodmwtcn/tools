function [potdens] = Z_fast_potdens(S,T)
% 1/30/2013 Parker MacCready
%
% Uses routines copied from the seawater routines to calculate
% potential density (for Pressure = 0).
%
% Faster becasue it works for any size of matrix
%
% Assumes T is potential temperature (deg C, P = 0)
% S = Salinity (same size as T)

%----------------------
% DEFINE CONSTANTS
%----------------------
a0 = 999.842594;
a1 =   6.793952e-2;
a2 =  -9.095290e-3;
a3 =   1.001685e-4;
a4 =  -1.120083e-6;
a5 =   6.536332e-9;

T68 = T * 1.00024;
smow = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68).*T68).*T68).*T68).*T68;

%     UNESCO 1983 eqn(13) p17.

b0 =  8.24493e-1;
b1 = -4.0899e-3;
b2 =  7.6438e-5;
b3 = -8.2467e-7;
b4 =  5.3875e-9;

c0 = -5.72466e-3;
c1 = +1.0227e-4;
c2 = -1.6546e-6;

d0 = 4.8314e-4;
potdens = smow + (b0 + (b1 + (b2 + (b3 + b4*T68).*T68).*T68).*T68).*S  ...
                   + (c0 + (c1 + c2*T68).*T68).*S.*sqrt(S) + d0*S.^2;


