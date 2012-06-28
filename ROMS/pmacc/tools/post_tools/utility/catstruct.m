function [z] = catstruct(x,y)
%CATSTRUCT  Concatenates the (matrix) elements of two structures.
%
%   edited by DAS, 4/20/2009 from JLab- now just concatenates vectors in y with those in x
%       - must have same fields!! 
%
%   Usage:  z=catstruct(x,y);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2003--2007 J.M. Lilly --- type 'help jlab_license' for details
%
% edited by N Banas jun 09 to deal gracefully with empty inputs

if isempty(x)
	z = y;
elseif isempty(y)
	z = x;
else
	  
	fx=fields(x);
	fy=fields(y);
	
	for k=1:length(fx)
	  if ~aresame(fx{k},fy{k})
		error('X and Y must have identical fields')
	  end
	end
	
	z=[];
	for i=1:length(fx);
	   vx=x.(fx{i});
	   vy=y.(fx{i});
	   %vx=getfield(x,fx{i});
	   %vy=getfield(y,fx{i});
	   N=max(size(vx,1),size(vy,1));
	 
	   if all(isreal(vx))  
		  vz=nan*zeros(N,size(vx,2)+size(vy,2));
	   else
		  vz=(nan+sqrt(-1)*nan)*zeros(N,size(vx,2)+size(vy,2));
	   end
	   %vz(1:size(vx,1),1:size(vx,2))=vx;
	   %vz(1:size(vy,1),[1:size(vy,2)]+size(vx,2))=vy;
	   %z=setfield(z,fx{i},vz);
	   %z.(fx{i})=vz;
	   z.(fx{i}) = [vx(:);vy(:)];
	end

end
