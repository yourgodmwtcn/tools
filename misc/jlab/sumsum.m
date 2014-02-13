function[b]=sumsum(x)
%SUMSUM(X)=SUM(X(ISFINITE(X)))
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2012 J.M. Lilly --- type 'help jlab_license' for details
b=sum(x(isfinite(x)));
