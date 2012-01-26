function t = isreal(x)
%ISREAL True for matrix that contains only real elements.
%       ISREAL(X) returns 1 if all elements in X have zero
%       imaginary part.

%       Clay M. Thompson 10-9-92
%       Copyright (c) 1992 by The MathWorks, Inc.
%       $Revision: 1.4 $  $Date: 1993/09/09 19:36:19 $

t = all(all(imag(x)==0));

