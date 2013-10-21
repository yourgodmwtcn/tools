function [Tide] =  PStide_at_seg( ISegment, StartDate, HDays, IP )
% [Tide] =  PStide_at_seg(ISegment,StartDate,HDays,IP)
% Computes hourly tidal height or current at specified segment
%
% Inputs:  
%	ISegment = segment number (cf. J.W. Lavelle et al., NOAA/PMEL
%		"A multiply-connected channel model of tides and tidal currents
%		in Puget Sound, WA and a comparison with updated observations",
%		NOAA Tech Memo. ERL PMEL-84, 1988.
%	StartDate = start date for time series, format must be 'dd-MON-yyyy' ;
%	HDays = number of days for time series;
%	IP = type of output:  1 = tidal heights, 2 = tidal currents;
% Output:
%	Tide = column vector of tidal heights(m) or currents(m/s),
%		starting at midnight (Pacific Standard Time) at the beginning of
%		date=StartDate and then hourly for HDays days.

%% This Matlab program file was adapted from the original FORTRAN code.
%% PStide_at_seg() is the computational part of SUBROUTINE TIDES();
%% Local functions Node() and Harmonic() replace respective SUBROUTINEs.
%% Original comments are retained (% ! ...), and much of the code was
%% kept similar to the FORTRAN for easier comparison.
%% Dave Winkel, APL, Jan 23, 2001

global PSfold % folder containing PS_TIDES programs and data files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  First, load in constituents for 22 tidal frequencies %%%
%    (Were in WASH.DAT)

load([PSfold 'PSTM_wash.mat'])
NAME = TIDEPARAM.name;
OMEGA = TIDEPARAM.omega; % degrees/hour
VU = TIDEPARAM.vu;
RELAMP = TIDEPARAM.relamp;
RELPHASE = TIDEPARAM.relphase;
MaxCon = length(OMEGA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Next, retrieve amplitude(factor) and phase(offset) corrections  %%%
%%%  for specific segment, at five dominant frequencies  

[O1,O1g,K1,K1g,N2,N2g,M2,M2g,S2,S2g] = Harmonic( IP,ISegment );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Now, compute Julian day,  and elasped days since Jan 1, 1985  %%%

tmp = datevec(StartDate);
IYear = tmp(1); % 4-digit year
elJ85 = datenum('01-Jan-1985');
elJyr = datenum(IYear,1,1);
eldt = datenum(StartDate);
IYd = eldt - elJyr + 1; % Julian day for year (Jan 1, IYear = 1)
IJd = eldt - elJ85 + 1; % Elapsed days since Jan 1, 1985 (Jan 1, 1985 = 1)

%%%%
% !  ******************************************************
% !     CALLS SUBROUTINE TO CALCULATE NODE FACTOR F FOR
% !     MIDDLE OF YEAR  %% (Actually, midpoint of time series ??)

[FNode] = Node ( IYd, IYear, HDays );

% !   *****************************************************
% !   CALCULATES MINOR CONSTITUENTS

G = NaN * OMEGA; H = G; % initialize

del1 = K1g - O1g;

for i=1:9
   H(i) = RELAMP(i) * O1;
   G(i) = K1g + RELPHASE(i) * del1;
end

H( 6 ) = RELAMP( 6 ) * K1;
H( 7 ) = K1;

delmn = M2g - N2g;
delsm = S2g - M2g;

for i=10:16
   H(i) = RELAMP(i) * M2;
   G(i) = M2g + RELPHASE(i) * delmn;
end

I = i+1; % To mimic Fortran do-loop on exit (I's here, or 11,15 ???)
% 11=MU 15=LAM 17=T2
G( 11 ) = S2g + RELPHASE( I ) * delsm;
G( 15 ) = S2g + RELPHASE( I ) * delsm;
H( 12 ) = N2;

for i=17:20
   H(i) = RELAMP(i) * S2;
   G(i) = S2g + RELPHASE(i) * delsm;
end
H( 23 ) = M2 + .688 * ( O1 + K1 ); % Mean height ??

H( 21 ) = 0.0;            % !SA
H( 22 ) = 0.0;            % !SSA
G( 21 ) = 0.0;            % !SAG
G( 22 ) = 0.0;            % !SSAG

% !     CONVERTING AMPLITUDES IN CM TO METERS
H = H / 100;

Mean = 0;
if IP==1
   Mean = H(23);
end
% !  CHANGE TO RADIANS     

Rad = pi / 180.0;

R = FNode .* H(1:22);
ETA = Rad * ( G - VU );
OMEGA = Rad * OMEGA;

% !  JULIAN DAYS CHANGED TO HOURS  %% (since 0000 1-1-1985 ??)
% !  CALCULATES NUMBER OF HOURS WANTING TO BE OBSERVED

jstart = 24 * ( IJd - 1 );
Length = 24 * HDays;

% !  CHANGE FROM GREENWHICH MEAN TIME TO PACIFIC STANDARD TIME

for j=1:Length
   timej = jstart + j - 1.0 + 8.0;       
   TideI = Mean;
   for i=1:MaxCon
      phase = OMEGA(i) * timej - ETA(i);
      TideI = TideI + R(i) * cos(phase);     
   end
   Tide(j,1) =  TideI;
end

return

% !   *********************************************************

function [FNode] =  Node( IYd, IYear, HDays )

EqCons = {'QQ';'Q1';'RHO';'O1';'M1';'P1';'K1';'J1';'OO';'NN';'MU'; ...
      'N2';'NU';'M2';'LAM';'L2';'T2';'S2';'R2';'K2';'SA';'SSA'};
for j=1:length(EqCons)
   FNode(j,1) = NaN;
end

% !  DM = THE DAY OF YEAR OF MIDPOINT OF SERIES
dm =  ( HDays / 2.0 ) + IYd;
% !  GRMS = GREENWICH HOUR AT THE BEGINNING OF THE SERIES
Grms = mod( HDays, 2 ) * 12.0 + 8.0;
% !   XX = THE CORRECTION FOR LEAP YEARS
xx = .25 * ( IYear - ( 1900 + 1.0 ) );

nprime = 259.18253333;
n = nprime - 19.32818576 * ( IYear - 1900 ) ...
   - .0529539336 * ( dm + xx ) - .0022064139 * Grms;  

r = 180.0 / pi;
x = n / r;
i = acos( .9136949 - (.0356926 *  cos(x)) );

% !  LG = LONGITUDE IN THE MOON'S ORBIT OF THE MOON'S INTERSECTION
% !      WITH THE CELESTIAL EQUATOR

lg = ( .206727 * sin(x) ) * ( 1.0 - .0194926 * cos(x) ); 
lg1 = .9979852 + .206727 * cos(x) - ( .0020148 * cos(2.0*x) );
lg2 = atan( lg / lg1 ) * r;

pprime = 334.32801944;

p = pprime + 40.66246584 * ( IYear - 1900.0 ) + ...
   .111404016 * ( dm + xx ) + .004641834 * Grms; 

pp = ( p - lg2 ) / r;

v = asin( .0897056 *  sin(x) / sin(i) );

% !  CALCULATES THE NODE FACTOR F FOR MIDDLE OF YEAR

O1 =  ( sin(i) * ( cos(.5*i) )^2 ) / .37988; 
K1 = 1 / ( ( (.8965 * sin(2.0*i)^2) + ...
   (.6001 * (sin(2.0*i) * cos(v)) ) + .1006 )^-.5 ); 

M2 = 1 / ( .91544 / (cos(.5*i )^4) );
K2 = 1 / ( (19.0444 * (sin(i)^4) + ...
   2.7702 * (sin(i)^2 * cos(2.0*v)) + .0981 )^-.5 );

J1 = 1 / ( .72137 / (sin(2.0*i)) );
OO = 1 / ( .016358 / (sin(i) * (sin(.5*i)^2)) );
L2 = M2 * ( ( 1.0 - 12.0 * (tan(.5*i )^2) * ...
   cos(2.0*pp) + 36.0 * (tan(.5*i)^4) )^.5 ); 
M1 = O1 * ( ( 2.310 + 1.435 * cos(2.0*pp) )^.5 ); 

Q1 = O1;
QQ = O1;
RHO = O1;

N2 = M2;
NN = M2;
LAM = M2;
MU = M2;
NU = M2;

S2 = 1.000;
R2 = 1.000;
T2 = 1.000;
P1 = 1.000;   
SA = 1.000;
SSA = 1.000;

for j=1:22
   eval( ['FNode(j) = ' EqCons{j} ';'] );
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [o1v,o1g,k1v,k1g,n2v,n2g,m2v,m2g,s2v,s2g] = Harmonic(ip,id);
% For segment # id, returns amplitude and phase adjustments at five frequencies.

% !        This subroutine accesses the direct access files containing the
% !        constituent values corresponding to a specified segment number and 
% !        returns the necessary constituent values to the WASHTIDE MODEL.

global PSfold

tide = ['o1';'k1';'n2';'m2';'s2'];
o1v=NaN; o1g=NaN; k1v=NaN; k1g=NaN; n2v=NaN; n2g=NaN; 
m2v=NaN; m2g=NaN; s2v=NaN; s2g=NaN; 

load([PSfold 'PSTM_seg_AmpPha.mat']) % replaces files=o1,k2,n2,m2,s2direct.dat
% holds structure AMP_PHA(1:5)

for ntide=1:5
   if ~strcmp( tide(ntide,1:2), AMP_PHA(ntide).name )
      error('PS, Harmonic(): Components not synchronized');
   end
   % vector index for segment should match segment number, but check anyway
   sn = AMP_PHA(ntide).segment;
   isn = find(sn==id);
   if isempty(isn)
      error(['PS, Harmonic(): Can''t find segment = ' num2str(id)]);
   elseif isn~=id
      warning(['PS, Harmonic(): index=' num2str(isn) ' for segment=' num2str(id)]);
   end
   % gather constituent values
   vel(ntide) = AMP_PHA(ntide).vel(isn);
   gphase(ntide) = AMP_PHA(ntide).gphase(isn);
   tamp(ntide) = AMP_PHA(ntide).tamp(isn);
   tphase(ntide) = AMP_PHA(ntide).tphase(isn);
end

if ip==1 % tidal heights
   o1v = tamp(1) * 100.0;
   o1g = tphase(1);
   k1v = tamp(2) * 100.0;
   k1g = tphase(2);
   n2v = tamp(3) * 100.0;
   n2g = tphase(3); 
   m2v = tamp(4) * 100.0;
   m2g = tphase(4);
   s2v = tamp(5) * 100.0;
   s2g = tphase(5);
else % area-averaged tidal flow
   o1v=vel(1)*100.0;
   o1g=gphase(1) + 180.0;
   k1v=vel(2)*100.0;
   k1g=gphase(2) + 180.0;
   n2v=vel(3)*100.0;
   n2g=gphase(3) + 180.0;
   m2v=vel(4)*100.0;
   m2g=gphase(4) + 180.0;
   s2v=vel(5)*100.0;
   s2g=gphase(5) + 180.0;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
