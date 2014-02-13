function[Ro]=ellrossby(lat,omega)
%ELLROSSBY  Ellipse Rossby number, for oceanographic applications.
%
%   RO=ELLROSSBY(LAT,OMEGA) returns the ellipse Rossby number RO for an 
%   ellipse with joint instantaeous frequency OMEGA, in radians per day, 
%   at latitude LAT.
%  
%   The ellipse Rossby number is defined as RO=2*OMEGA/F where F is the
%   Coriolis frequency at latitude LAT.  
%
%   LAT and OMEGA may be cell arrays of numeric arrays, both having the
%   same size.  RO will then be a similar cell array.
%
%   Usage: ro=ellrossby(lat,omega);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    ellrossby_test,return
end

if ~iscell(lat)
   Ro=ellrossby_one(lat,omega);
   Ro(isinf(lat.*omega))=inf;
else
    for i=1:length(lat)
          Ro{i,1}=ellrossby_one(lat{i},omega{i});
          Ro{i,1}(isinf(lat{i}.*omega{i}))=inf;
    end
end


function[Ro]=ellrossby_one(lat,omega)
fcor=2*pi*corfreq(lat)*24;  %Coriolis frequency in cycles per day at center 
Ro=frac(2*omega,fcor);       %Rossby number under solid-body assumptio
    
    
function[]=ellrossby_test
 
%reporttest('ELLROSSBY',aresame())
