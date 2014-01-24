function [dist, distFrom]=trackDist(varargin)
% [dist, distFrom]=trackDist(x,y,xTrack,yTrack);
%                         ...x,y,track);
%                         ...EX,track);
%                         ...EX,xTrack,yTrack);
%
% returns distance along and distance from a trackline.
% if x or y has nans, dist and distFrom will have nans in those locations


if nargin==2
    Ex=varargin{1};
    x=Ex.x;
    y=Ex.y;
    track=varargin{2};
    xTrack=track.x;
    yTrack=track.y;
    
elseif nargin==4
    x=varargin{1};
    y=varargin{2};
     xTrack=varargin{3};
    yTrack=varargin{4};
    
elseif nargin==3
    if isstruct(varargin{1})==1
        Ex=varargin{1};
    x=Ex.x;
    y=Ex.y;
    xTrack=varargin{2};
    yTrack=varargin{3};
    else
    x=varargin{1};
    y=varargin{2};
    track=varargin{3};
    xTrack=track.x;
    yTrack=track.y;
    end
    

end

%% check that track points do not have nans or bad points, if so, delete points and return a warning
badtrack=find(isnan(xTrack)==1 |isnan(yTrack)==1);
if isempty(badtrack)==0
xTrack(badtrack)=[];
yTrack(badtrack)=[];
warning('Section has nan data points, these section points will not be used in determining distance along track')
end

%% change to x y coordinate space

lat=y*111.325;
lon=cosd(mean(yTrack))*111.325.*x;



tlat=yTrack*111.325;
tlon=cosd(mean(yTrack))*111.325.*xTrack;
%% special case where only 1 track point is given, return dist=0 and distFrom=?
if length(tlon)==1

distFrom=sqrt((tlon-lon).^2+(tlat-lat).^2);
dist=zeros(length(lon),1);
return

end

%%

%first find distance between segments if not there
if exist('track','var')==1  
    if  isfield(track, 'dist')
    Tdist=track.dist;
    elseif  isfield(track, 'along_dist')
    Tdist=track.along_dist;
    else
    TrackDistance=zeros(length(tlon),1);
        for k=1:length(tlon)-1;
        TrackDistance(k+1)=sqrt((tlon(k)-tlon(k+1)).^2+(tlat(k)-tlat(k+1)).^2);
        end
        Tdist=cumsum(TrackDistance); 
    end
    
else   
    TrackDistance=zeros(length(tlon),1);
        for k=1:length(tlon)-1;
        TrackDistance(k+1)=sqrt((tlon(k)-tlon(k+1)).^2+(tlat(k)-tlat(k+1)).^2);
        end
        Tdist=cumsum(TrackDistance);

   
end
%%

Tdist_point=nan(length(lon),1);
MinDistance=nan(length(lon),1);
    %put 99999 in places where lon or lat is currently nan and replace at the
    %end with nans

    nanlonlat=find(isnan(lon)==1 | isnan(lat)==1);
    Tdist_point(nanlonlat)=9999;
    MinDistance(nanlonlat)=9999;
q=1;
while q>0

    

    if MinDistance(q)~=9999
    distToPoint=sqrt((tlon-lon(q)).^2+(tlat-lat(q)).^2);
    [D,I]=min(distToPoint);
 

    pl=find(lon==lon(q) & lat==lat(q));
    
    MinDistance(pl,1)=D(1) ;
    
    Tdist_point(pl)=Tdist(I(1));
    end
    
     nanfind=find(isnan(MinDistance(:,1))==1,1);
    
          if isempty(nanfind)==0
            q=nanfind(1);
          else q=0;
          end
          
end
Tdist_point(nanlonlat)=nan;
MinDistance(nanlonlat)=nan;
dist=Tdist_point;
distFrom=MinDistance(:,1);

%%
