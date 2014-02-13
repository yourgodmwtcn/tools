function[varargout]=cell2col(varargin)
%CELL2COL  Converts cell arrays of column vectors into 'column-appended' data.
%
%   COL=CELL2COL(X) concatenates X, a cell array of N column vectors, into a 
%   single column vector COL with N blocks of data, each ending with a NAN.
%
%   To mark missing data, any pre-exiting NANs in X will be replaced with 
%   INFs, and any complex NANs will be replaced with complex INFs.
%
%   CELL2COL(X) where X is a cell array of numeric arrays of any dimension
%   also works.  In this case each array will be convered into a column 
%   vector via X{1}=X{1}(:) and so forth prior to concatenation.
%   __________________________________________________________________
%
%   Multiple input /output arguments
%
%   [C1,C2,...,CN]=CELL2COL(X1,X2,...,XN) also works.  
%   __________________________________________________________________
% 
%   Variable overwriting
%
%   CELL2COL(X1,X2,...,XN) with no output arguments overwrites the input
%   variables.
%   __________________________________________________________________
%
%   Invertibility
%  
%   CELL2COL is inverted by COL2CELL, provided 
%
%       (i)  The cell arrays XN are all column vectors of the same size and
%       (ii) The first input argument X1 contains no NANs.
%   __________________________________________________________________
%
%   See also COL2CELL, COL2MAT, MAT2COL, COLBREAKS, VCELLCAT.
%
%   Usage: col=cell2col(x);
%          [c1,c2,...,cN]=cell2col(x1,x2,...,xN);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2012 J.M. Lilly --- type 'help jlab_license' for details


for i=1:nargin
    ray=varargin{i}(:);  
    %if 1
    for j=1:length(ray)
        ray{j}=[ray{j};-inf];
    end
    ray=cell2mat(ray);
    
    ray=vswap(ray,nan,inf);
    ray=vswap(ray,nan+sqrt(-1)*nan,inf+sqrt(-1)*inf);
    ray=vswap(ray,-inf,nan);

    if i==1
        ray1=ray;
    elseif i~=1
        ray(isnan(ray1))=nan;
        ray(isinf(ray1))=nan;
    end
        
    varargout{i}=ray;
%   Former version, my way, not as fast
%     elseif 0
%     for j=1:length(ray)
%         temp=ray{j}(:);
%         %Replace missing data i
%         temp=vswap(temp,nan,inf);
%         temp=vswap(temp,nan+sqrt(-1)*nan,inf+sqrt(-1)*inf);
%         ray{j}=[temp;nan];
%         %Make sure pattern of NANs matches that in the first argument, the key
%         if i~=1
%             ray{j}(isnan(ray1{j}))=nan;
%             ray{j}(isinf(ray1{j}))=nan;
%         end 
%     end
%     if i==1
%         ray1=ray;
%     end
%     varargout{i}=vcellcat(ray);
%     end
end   
eval(to_overwrite(nargin));

function[]=cell2col_test
%Tests for cell2col are contained in col2cell