function[h]=cellplot(varargin)
%CELLPLOT  Plot elements of a cell array of numeric arrays.
%
%   CELLPLOT(X) where X is a cell array containing N arrays,
%
%       X{1}=X1, X{2}=X2,..., X{N}=XN
% 
%   is the same as 
%
%      PLOT(X1), PLOT(X2),..., PLOT(XN).
%
%   However, CELLPLOT is vastly quicker than the obvious loop and also has
%   some handy features, as described below.
%
%   CELLPLOT(X,Y), where X and Y are both cell arrays of N arrays, performs
%
%      PLOT(X1,Y1), PLOT(X1,Y1),..., PLOT(XN,YN).
%
%   H=CELLPLOT(...) returns the handle to the plots.  
%   __________________________________________________________________
%   
%   Additional options
%
%   CELLPLOT(...,STY) uses the linestyle specified by STY.  STY is a
%   string following the format in LINESTYLE, e.g. STY='2b g r--'. 
%
%   CELLPLOT(...,'M_MAP') will work with Rich Pawlowicz's M_MAP 
%   package by calling M_PLOT.  
%
%   CELLPLOT(...,INDEX) only plots the elements of the cell array 
%   indicated by index.
%
%   The string arguments and INDEX can be combined provided they are 
%   after the one or two input cell array variables.
%   __________________________________________________________________
%
%   CELLPLOT on the sphere
%    
%   When plotting data on the sphere, there is an annoying wrap-around
%   effect when the data crosses the dateline at longitude 180.  This
%   leads to lines extending all the way across the plot.
%
%   CELLPLOT(LONO,LON,LAT) where LONO is the value of a cutoff longitude
%   uses LONO as the right-hand-edge of the plot, and LONO-360 as the 
%   left-hand-edge, with no wraparound effects.  
%
%   This works together also with the 'm_map' option.  In this case LONO 
%   should be between -180 and 180, and one must set the maximum longitude
%   in M_MAP to LONO. For a global plot, this means one will call M_PROJ as
% 
%         M_PROJ(...,'longitudes',[LONO-360 LONO],...).
%
%   If you run CELLPLOT with M_MAP and data doesn't show up, make sure 
%   the longitudes have been set correctly.
%   __________________________________________________________________
%
%   See also LINESTYLE.
%
%   Usage: cellplot(x)
%          cellplot(x,y)
%          h=cellplot(x,y,index,'m_map');
%          cellplot(lono,lat,lon,index);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2012 J.M. Lilly --- type 'help jlab_license' for details
 
linestr=[];
str='matlab';

lono=[];
if ~iscell(varargin{1})
    if length(varargin{1})==1
        lono=varargin{1};
        varargin=varargin(2:end);
    end
end

while ischar(varargin{end})
    if ~isempty(strfind(varargin{end},'mat'))||~isempty(strfind(varargin{end},'m_m'))
        str=varargin{end};
    else
        linestr=varargin{end};
    end
    varargin=varargin(1:end-1);
end

if  ~iscell(varargin{end}) 
    ii=varargin{end};
    varargin=varargin(1:end-1);
else
    ii=inf;
end
holdstate=ishold;
storestate=get(gcf,'BackingStore');

if ~holdstate
    cla
end

if length(varargin)==0
    error('X is not a cell array.')
end
if ~iscell(varargin{1})
    error('X is not a cell array.')
end

x=varargin{1};
if length(varargin)==1
     y=[];
else
     y=varargin{2}; 
end

if isempty(ii)
    x=x(ii);
else
    if isfinite(ii)
        x=x(ii);
        if ~isempty(y)
            y=y(ii);
        end
    end
end

if ~isempty(x)
    if isempty(y)
        cell2col(x);
        col2mat(x);
        h=plot(x);
    else
        cell2col(x,y);
        
        if ~isempty(lono)
             %maxmax(frac(360,2*pi)*angle(rot(jdeg2rad(x)-jdeg2rad(lono))))
             %x=frac(360,2*pi)*unwrap(angle(rot(frac(2*pi,360)*(x-lono))))+lono;
             x=deg360(frac(360,2*pi)*unwrap(angle(rot(frac(2*pi,360)*(x-lono)))))+lono;
             for i=1:size(x,2)
                 while maxmax(x(:,i))>lono
                    x(:,i)=x(:,i)-360;
                 end
             end
             %x=deg360(x-deg360(lono))+deg360(lono)-360;
             %booljump=abs(x-vshift(x,-1,1))>90;
             booljump=abs(x(2:end)-x(1:end-1))>90;
             x(booljump)=nan;
             y(booljump)=nan;
             col2mat(x,y);  
        else
             %x(end+1,:)=nan;
             %y(end+1,:)=nan;
             col2mat(x,y);
        end

        if strcmp(str(1:3),'mat')
              h=plot(x,y);hold on
        elseif  strcmp(str(1:3),'m_m')
              h=m_plot(x,y);hold on
        end
    end
end

if ~isempty(linestr)
    linestyle(h,linestr);
end

if nargout==0
    clear h
end


set(gcf,'BackingStore',storestate)

if ~holdstate
    hold off
end

