function[region]=latlonregion(lat,lon,delta)
%LATLONREGION
%
%   LATLONREGION
%
%   'latlonregion --t' runs a test.
%
%   Usage: []=latlonregion();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(lat, '--t')
    latlonregion_test,return
end

if nargin==2
    delta=0;
end

if ~iscell(lat)
    region=latlonregion_one(lat,lon,delta);
else
    [latm,lonm]=latlonregion_mid(lat,lon,delta);
    for i=1:length(lat)
        regioni(i,:)=latlonregion_one(lat{i},deg180(lon{i}-lonm),delta);
    end
    regioni;
    region(1)=minmin(regioni(:,1))+lonm;
    region(2)=maxmax(regioni(:,2))+lonm;
    region(3)=minmin(regioni(:,3));
    region(4)=maxmax(regioni(:,4));
end

        
function[region]=latlonregion_one(lat,lon,delta)


minlat=max(minmin(lat)-delta,-90);
maxlat=min(maxmax(lat)+delta,90);

lonunwrap=frac(360,2*pi)*unwrap(angle(rot(frac(2*pi,360)*lon)));
%figure,plot(lonunwrap,lat)
dlon=maxmax(lonunwrap)-minmin(lonunwrap)+2*delta;
minlon=deg180(minmin(lonunwrap))-delta;
maxlon=minlon+dlon;

region=[minlon maxlon minlat maxlat];

%First estblish the "center" based on mean value

function[latm,lonm]=latlonregion_mid(lat,lon,delta)

[x,y,z]=latlon2xyz(cell2col(lat),cell2col(lon));
[mx,my,mz]=vmean(x(:),y(:),z(:),1);
alpha=sqrt(mx.^2+my.^2+mz.^2)./radearth;
mx=mx./alpha;my=my./alpha;mz=mz./alpha;
[latm,lonm]=xyz2latlon(mx,my,mz);


% lonm
% figure,plot(angle(rot(jdeg2rad(lon-lonm))))
% frac(360,2*pi)*minmin(angle(rot(jdeg2rad(lon-lonm))))

%figure,plot(frac(360,2*pi)*angle(rot(frac(2*pi,360)*(lon-lonm)))+lonm)


%minlon=minmin(frac(360,2*pi)*angle(rot(frac(2*pi,360)*(lon-lonm)))+lonm);
%maxlon=maxmax(frac(360,2*pi)*angle(rot(frac(2*pi,360)*(lon-lonm)))+lonm);

%minlon=lonm+frac(360,2*pi)*minmin(angle(rot(jdeg2rad(lon-lonm))))-delta;
%maxlon=lonm+frac(360,2*pi)*maxmax(angle(rot(jdeg2rad(lon-lonm))))+delta;
 
% if minlon<-180
%     minlon=minlon+360;
%     maxlon=maxlon+360;
% end



function[]=latlonregion_test



%reporttest('LATLONREGION',aresame())
