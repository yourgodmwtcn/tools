%/************************************************************************************************    
%Analysis for Amy Bower's expath floats
str='noprint';
load EXPfloatmaster;
lat=Y;
lon=X;

num=Jdates'-2440000+datenum(1968,5,23,0,0,0); %Jdate 244000 = 0000 hours, May 23, 1968, see "gregorian"
cv=latlon2uv(num,lat,lon);
id=flts;
t=T;
p=P;

%581,664,680;
%jj=find(id==581);

lon=fillbad(lon);
lat=fillbad(lat);

matsave expfloat num lat lon id t p cv 

use expfloat

jj=find(id==664);
%jj=find(id==581);
%jj=find(id==680);
vindex(lat,lon,id,t,p,cv,jj,2);
ii=min(find(~isnan(lat))):max(find(~isnan(lat)));
vindex(num,lat,lon,id,t,p,cv,ii,1);

%plot(lon,lat)

lato=vmean(lat,1);
lono=vmean(lon,1);

cx=latlon2xy(lat,lon,lato,lono);
fcor=2*pi*corfreq(lat)*24; 
 
dt=1;
ga=3;
be=7;

fs=morsespace(ga,be,{0.05,pi/2},2*pi/20,8);
[wx,wy]=wavetrans(real(cx),imag(cx),{1,ga,be,fs,'bandpass'},'mirror');

%Ridges
[ir,jr,xr,fr,br,cr]=ridgewalk(1,wx,wy,fs,{morseprops(ga,be),2,'amp'});  
[d1r,d2r]=mom2dev(xr,fr,br,cr,2);
omegar=powermean(fr,xr,2);
xerrr=frac(1,2)*squared(frac(sqrt(ga*be),vrep(omegar,2,2))).*d2r;

[xhat,yhat,omega,xerr,yerr]=ridgemap(length(lat),xr(:,1),xr(:,2),omegar,xerrr(:,1),xerrr(:,2),ir,'collapse');



%Ellipse properties
[kappa,lambda,theta,phi]=ellparams(xhat,yhat);
R=ellrad(kappa,lambda);
V=ellvel(dt*24*3600,kappa,lambda,theta,phi,1e5); %Centimeters per second
Ro=frac(2*omega,fcor);  %Rossby number under solid-body assumption
xi=sign(lambda).*sqrt(1-lambda.^2);
xhat=real(xhat)+sqrt(-1)*real(yhat);

%Error properties
[kaperr,lamerr,theerr,phierr]=ellparams(xerr,yerr);

%Residual
[lathat,lonhat]=xy2latlon(real(xhat),imag(xhat),lato,lono);
latres=lat-(vswap(lathat-lato,nan,0));
lonres=lon-(vswap(lonhat-lono,nan,0));


%make struct.id664 num lat lon latres lonres xhat kappa xi theta phi omega R V Ro kaperr lamerr theerr phierr
%make struct.id581 num lat lon latres lonres xhat kappa xi theta phi omega R V Ro kaperr lamerr theerr phierr
make struct.id680 num lat lon latres lonres xhat kappa xi theta phi omega R V Ro kaperr lamerr theerr phierr

%\************************************************************************************************    


%use struct.id664
use struct.id680


%/*************************************************    
%Figure....
use expfloat

%jj=find(id==664);
%jj=find(id==581);
jj=find(id==680);
vindex(lat,lon,id,t,p,cv,jj,2);
ii=min(find(~isnan(lat))):max(find(~isnan(lat)));
vindex(num,lat,lon,id,t,p,cv,ii,1);

%plot(lon,lat)

lato=vmean(lat,1);
lono=vmean(lon,1);

cx=latlon2xy(lat,lon,lato,lono);
fcor=2*pi*corfreq(lat)*24; 
 
dt=1;
ga=3;
be=7;

fs=morsespace(ga,be,{0.05,pi/2},2*pi/20,8);
[wx,wy]=wavetrans(real(cx),imag(cx),{1,ga,be,fs,'bandpass'},'mirror');

%Ridges
[ir,jr,xr,fr,br,cr]=ridgewalk(1,wx,wy,fs,{morseprops(ga,be),2,'amp'});  
[d1r,d2r]=mom2dev(xr,fr,br,cr,2);

%hold on,plot(num-num(1),2*pi./omega,'r')
%\*************************************************    


  
%/*************************************************    
figure
ci=(0:2:100);
numo=num(1);
[h,hl]=wavespecplot(num-numo,cv,2*pi./fs,sqrt(abs(wp).^2+abs(wn).^2),1,ci);


linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num{ii}-numo,2*pi./fbar{ii},'w','linewidth',4)
plot(num{ii}-numo,2*pi./fbar{ii},'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')


orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])
if strcmp(str,'print')
    print -deps ebasn-example.eps
end
%\********************



struct=ellipsextract(num,lat,lon,4.5,2,1/2,1/64);


%/***********************************************************
%Fake eddy plots
num=[-50:1/24:50]';

%10 km radius eddy that starts at time num=0
x=[10+num/5].*rot(-2*pi*num/5);
x(num<0)=0;

Rtrue=abs(x);
omtrue=2*pi/5+0*x;
Vtrue=-Rtrue.*omtrue*frac(1000*100,24*3600);
xitrue=1+0*x;xitrue(num<0)=0;

%Convert to lat/lon at 50 degrees north
[lat,lon]=xy2latlon(real(x),imag(x),50,0);

struct=eddyextract(num,lat,lon,3,1,1/2,1/64);

use struct
figure,
subplot(4,1,1),
uvplot(num(num<5),xhat(num<5)),hold on
uvplot(num(num>45),xhat(num>45)),
uvplot(num(num<45&num>5),xhat(num<45&num>5)),
uvplot(num,x),uvplot(num,x),
ylabel('Displacement (km)'),ylim([-12 12]*2),hold on
linestyle 2G-- 2G-- 2G-- 2G-- 2G 2G 2w 2w k k--
title('Analysis of a Synthetic "Eddy"')
subplot(4,1,2),
plot(num(num<5),[V(num<5) Vtrue(num<5)]),hold on
plot(num(num>45),[V(num>45) Vtrue(num>45)]),
plot(num(num<45&num>5),[V(num<45&num>5) Vtrue(num<45&num>5)]),
plot(num(num<45&num>5),[V(num<45&num>5) Vtrue(num<45&num>5)]),
ylabel('V (cm/s)'),ylim([-32 0])
subplot(4,1,3),
plot(num(num<5),[R(num<5) Rtrue(num<5)]),hold on
plot(num(num>45),[R(num>45) Rtrue(num>45)]),
plot(num(num<45&num>5),[R(num<45&num>5) Rtrue(num<45&num>5)]),
plot(num(num<45&num>5),[R(num<45&num>5) Rtrue(num<45&num>5)]),
ylabel('R (km)'),ylim([0 22]),hold on
subplot(4,1,4)
plot(num(num<5),2*pi./[omega(num<5) omtrue(num<5)]),hold on
plot(num(num>45),2*pi./[omega(num>45) omtrue(num>45)]),
plot(num(num<45&num>5),2*pi./[omega(num<45&num>5) omtrue(num<45&num>5)]),
plot(num(num<45&num>5),2*pi./[omega(num<45&num>5) omtrue(num<45&num>5)]),
ylabel('Period 2 \pi/|\omega| (days)'),ylim([0 9.5])
%subplot(5,1,5),plot(num,[xi xitrue]),ylabel('Circularity Xi'),ylim([0 1.1]),xlabel('Days from Start of Eddy')
xlabel('Time (days)')

for i=1:4
    subplot(4,1,i),if i>1, linestyle 2G-- k 2G-- k 2G 2w 2G k, end
    xlim([-9 50])%,vlines(0,'G')
    vlines([5 45],'G:')
end
letterlabels(4)

packrows(4,1)
fontsize 12 10 10 10
orient tall
set(gcf,'paperposition',[1 1 3.5 8])
print -depsc fakeeddy.eps


figure
plot(R(num<5),V(num<5)),hold on
plot(R(num>45),V(num>45))
plot(R(num<45&num>5),V(num<45&num>5))
plot(R(num<45&num>5),V(num<45&num>5))
plot(Rtrue(num>0),Vtrue(num>0)),hold on
plot(Rtrue(num>0),Vtrue(num>0)),hold on
linestyle 3G-- 3G-- 4G 4G 3w 2k
title('Radius / Velocity Plot for Synthetic "Eddy"')
xlabel('Geometric Mean Radius (km)')
ylabel('Geometric Mean Velocity (cm/s)')
axis([0 21 -31 0])

fontsize 10 8 8 8 
orient portrait
set(gcf,'paperposition',[1 1 3.5 3.25])
print -depsc fakeeddy_rv.eps

%\***********************************************************





