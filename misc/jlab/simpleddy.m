function[v,psi,zeta,psif]=simpleddy(varargin)
%SIMPLEDDY  Streamfunction, velocity, and vorticity for various eddy profiles. 
%
%   [V,PSI,ZETA]=SIMPLEDDY(R,RO,VO) returns the streamfunction PSI,
%   azimuthal velocity V, and relative vorticity ZETA at radial locations R
%   for the Rankine vortex.
%
%   The eddy core radius is RO, in km, and at this radius the eddy 
%   currents achieve their maximum value of |VO| cm/s.  The sign of VO
%   gives the rotation sense of the eddy, positive for counterclockwise.
%
%   PSI has units of km^2/s, V has units of cm/s, and ZETA of 1/s. 
%
%   R should be a column vector.  RO and VO may either be scalars or arrays
%   of the same length N. The output arrays will have size LENGTH(R) x N.
%   _______________________________________________________________________
%
%   Eddy models
%
%   SIMPLEDDY can also return profiles associated with various other types
%   of model eddies.  The following forms are available:
%
%         SIMPLEDDY(R,RO,VO,'rankine')          --Rankine vortex  [default]
%         SIMPLEDDY(R,RO,VO,'gaussian')         --Gaussian vortex
%
%   Futher types of eddy models will be available in a future release. 
%   _______________________________________________________________________
% 
%   Eddy 'slices'
%
%   SIMPLEDDY can also be used to return an eddy 'slice', that is, the
%   streamfunction and vector-valued currents due to a moving eddy.  
%  
%   [CV,PSI,ZETA]=SIMPLEDDY(CX,RO,VO,...) where CX=X+iY is the complex-
%   valued position of the eddy center, returns the complex-valued eddy 
%   currents CV=U+iV observed at the origin.  
%   _______________________________________________________________________
%
%   'simpleddy --f' makes some sample figures.
%
%   Usage: [v,psi,zeta]=simpleddy(r,RO,VO);
%          [v,psi,zeta]=simpleddy(r,RO,VO,'gaussian');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
%         SIMPLEDDY(R,RO,VO,'matern',ALPHA)     --Matern vortex
%         SIMPLEDDY(R,RO,VO,'morse',GAMMA)      --Morse vortex
%
%   The Matern and Morse wavelet families are described by additional 
%   shape parameters, termed ALPHA and GAMMA respectively.  In order for
%   the profiles to be well defined, we must have ALPHA>1/2 and GAMMA>1.
%
%   RO, VO, and the shape parameter ALPHA or GAMMA may then each be either 
%   scalars or arrays of the same length. 

%   See also MATERNFUN, EDDYCURVE, EDDYFIT, EDDYGUESS.


%   XXX what if X has more than one column.  
%          [v,psi,zeta]=simpleddy(r,RO,VO,'matern',ALPHA);


if strcmp(varargin{1}, '--f')
    gaussian_fig
    %materneddy_figures
    %morseddy_figures
    return
end

x=varargin{1};
if size(x,1)==1
    disp('Transposing the first argument to SIMPLEDDY to make it a column vector.')
    x=conj(x');
end

ro=varargin{2};
vo=varargin{3};
if nargin==3
    str='ran';
else
    str=lower(varargin{4});
end
if strcmp(str(1:3),'mat')
    [v,psi,zeta]=materneddy(varargin{5},ro,vo,abs(x));
elseif strcmp(str(1:3),'gau')
    [v,psi,zeta]=gaussianeddy(ro,vo,abs(x));
elseif strcmp(str(1:3),'ran')
    [v,psi,zeta]=rankineddy(ro,vo,abs(x));
elseif strcmp(str(1:3),'mor')
    [v,psi,zeta]=morseddy(varargin{5},ro,vo,abs(x));
end


if size(x,2)==1
    x=vrep(x,size(v,2),2);
end
if ~isreal(x)
    v=-sqrt(-1)*frac(x,abs(x)).*v;
end

function[v,psi,zeta]=rankineddy(ro,vo,r)
%RANKINEDDY Velocity and streamfunction for a Rankine vortex.
%
%   RANKINEDDY computes the velocity and streamfunction observed at
%   the origin due to a Rankine vortex located at a specified point,
%   as described in Lilly and Rhines (2002).
%
%   [v,psi]=rankineddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'rankineddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    


vo=vo/100/1000;
iO=find(r>ro);


v=vo.*r./ro;
vO=vo.*ro./r;
if ~isempty(iO)
 %size(r),size(v), size(ro),size(iO), size(vo)
    v(iO)=vO(iO);
end
v=v*100*1000;

psi=(1/2).*(r./ro).^2;
psiO=log(r./ro)+1/2;
if ~isempty(iO)
    psi(iO)=psiO(iO);
end
psi=psi.*vo.*ro;

zeta=2*vo./ro+0*r;
if ~isempty(iO)
    zeta(iO)=0;
end

function[]=rankineddy_fig
ro=10;
vo=10;
eta=sqrt(-1)*5+(-50:.1:50)';
[v,psi]=rankineddy(eta,ro,vo);

figure,
subplot(121),
uvplot(real(eta),v),
title('Cyclone sliced south of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Rankine eddy')


function[v,psi,zeta]=gaussianeddy(ro,vo,r)
%GAUSSIANEDDY Velocity and streamfunction for a Gaussian vortex.
%
%   GAUSSIANEDDY computes the velocity and streamfunction observed at
%   the origin due to an eddy having a Gaussian streamfunction profile
%   located at a specified point.
%
%   [v,psi]=gaussianeddy(eta,ro,vo);
%
%     Input
%	eta: complex-valued position of eddy (km)
%	ro:  size of eddy (km)
%	vo:  velocity at edge, negative for anticyclone (cm/s) 
%
%     Output
%	v,psi:	azimuthal velocity and streamfunction at origin
%
%   'gaussianeddy --f' generates a sample figure
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2012 J.M. Lilly --- type 'help jlab_license' for details    

vo=vo./100./1000;

psi=-1.*vo.*ro.*exp(1/2.*(1-r.^2./ro.^2));
v=100*1000*r.*(vo./ro).*exp(1/2.*(1-r.^2./ro.^2));
zeta=(vo./ro).*(2-r.^2./ro.^2).*exp(1/2.*(1-r.^2./ro.^2));


function[]=gaussian_fig
ro=10;
vo=10;
eta=sqrt(-1)*5+(-50:.1:50)';

[v,psi]=gaussianeddy(eta,ro,vo);

figure,
subplot(121),
uvplot(real(eta),v),
title('Cyclone sliced south of center, halfway to edge')
subplot(122)
polar(angle(v),abs(v))
title('Hodograph of Gaussian eddy')


function[v,psi,zeta]=materneddy(alpha,R,V,r)
%MATERNEDDY
%
%   MATERNEDDY
%
%   Usage: []=materneddy();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011--2012 J.M. Lilly --- type 'help jlab_license' for details
 
V=V./100./1000;


Nalpha=length(alpha(:))';
NR=length(R(:))';
NV=length(V(:))';
N=max([Nalpha NR NV]);

alpha(alpha==1/2)=1/2+1/1000;

if N>1
    if Nalpha~=N
        if Nalpha==1
            alpha=alpha+zeros(1,N);
        else
            error('ALPHA must be scalar or an array of the same length as RO and VO.')
        end
    end
    if NR~=N
        if NR==1
            R=R+zeros(1,N);
        else
            error('RO must be scalar or an array of the same length as ALPHA and VO.')
        end
    end
    if NV~=N
        if NV==1
            V=V+zeros(1,N);
        else
            error('VO must be scalar or an array of the same length as RO and ALPHA.')
        end
    end
end

%Halpha=gamma(alpha).*(2.^(alpha-1));  
%Don't need this because I only need the ratio H/Halpha
c=-frac(maternrad(alpha).^2,R).*maternfun(alpha-1,maternrad(alpha));
Hnorm=frac(V,c);
Rnorm=R./maternrad(alpha);


if ~aresame(size(alpha),size(r))
    if size(r,2)==1&&(length(alpha)~=1)
        r=vrep(r,length(alpha),2); 
    end
    if size(r,1)~=1&&(length(alpha)~=1)
        Rnorm=vrep(Rnorm,size(r,1),1);
        Hnorm=vrep(Hnorm,size(r,1),1);   
    end
end
v=-100*1000*Hnorm.*r.*(1./Rnorm.^2).*maternfun(alpha-1,r./Rnorm);

if nargout >=2
    psi=Hnorm.*maternfun(alpha,r./Rnorm);
end
if nargout ==3
    zeta= Hnorm.*(1./Rnorm.^2).*((r./Rnorm).^2.*maternfun(alpha-2,r./Rnorm)-2*maternfun(alpha-1,r./Rnorm));
end



function[v,psi,zeta]=morseddy(gamma,R,V,r)
%MORSEDDY
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012 J.M. Lilly --- type 'help jlab_license' for details
 
V=V./100./1000;


Ngamma=length(gamma(:))';
NR=length(R(:))';
NV=length(V(:))';
N=max([Ngamma NR NV]);

if N>1
    if Ngamma~=N
        if Ngamma==1
            gamma=gamma+zeros(1,N);
        else
            error('GAMMA must be scalar or an array of the same length as RO and VO.')
        end
    end
    if NR~=N
        if NR==1
            R=R+zeros(1,N);
        else
            error('RO must be scalar or an array of the same length as GAMMA and VO.')
        end
    end
    if NV~=N
        if NV==1
            V=V+zeros(1,N);
        else
            error('VO must be scalar or an array of the same length as RO and GAMMA.')
        end
    end
end

Hnorm=-V.*R.*frac(1,gamma-1).*exp(frac(gamma-1,gamma));
Rnorm=R;


if ~aresame(size(gamma),size(r))
    if size(r,2)==1&&(length(gamma)~=1)
        r=vrep(r,length(gamma),2); 
    end
    if size(r,1)~=1&&(length(gamma)~=1)
        Rnorm=vrep(Rnorm,size(r,1),1);
        Hnorm=vrep(Hnorm,size(r,1),1);   
        gamma=vrep(gamma,size(r,1),1);   
    end
end

psi=Hnorm.*exp(-frac(gamma-1,gamma).*(r./Rnorm).^gamma);
v=-100*1000*frac(gamma-1,Rnorm).*(r./Rnorm).^(gamma-1).*psi;

if nargout ==3
    zeta= squared(frac(gamma-1,Rnorm)).*(r./Rnorm).^(gamma-2).*psi...
    .*((r./Rnorm).^gamma-frac(gamma,gamma-1));
end

function[]=materneddy_test
 
r=[0:.1:50]';
alpha=[1/2:1/2:10];
psi=materneddy(alpha,10,15,r);

function[]=materneddy_figures
dr=0.1;
r=[0:dr:500]';
%alpha=[1/2:1/4:10];
%alpha=logspace(log10(1.05),log10(25),10);
alpha=logspace(log10(3/4),log10(50),10);
alpha=([1.2:.2:2 3:2:11]);
[v,psi,zeta]=simpleddy(r,10,15,'matern',alpha);
[vo,psio,zetao]=simpleddy(r,10,15,'gaussian');
%[vo,psio,zetao]=simpleddy(r,10,15,'rankine');

psi=[psio psi];
v=[vo v];
zeta=[zetao zeta];

c=frac(maternfun(alpha-1,maternrad(alpha)),maternfun(alpha-1,0));
c=[1 c.*sqrt(exp(1))];
for i=1:length(c)
    psi(:,i)=psi(:,i).*c(i);
    v(:,i)=v(:,i).*c(i);
    zeta(:,i)=zeta(:,i).*c(i);
end

psitwosides=[flipud(psi);psi(2:end,:)];
Psi=abs(fft(psitwosides));
%[mu,sigma,skew,kurt]=pdfprops(1:length(Psi),psitwosides(:,[2:end 1]));
%[muf,sigmaf,skewf,kurtf]=pdfprops(1:length(Psi),fftshift(Psi(:,[2:end 1])));

f=fourier(pi,length(Psi));
Psi=Psi(1:length(r),:);

figure
subplot(2,2,1),h=plot(r,psi);vlines(10,'k:'), xlim([0 50]),
h1=h(1);h=h(2:end);title('Streamfunction'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,2),h=plot(r,v);vlines(10,'k:'),ylim([0 15.5]), xlim([0 50])
h1=h(1);h=h(2:end);title('Velocity'),xlabel('Distance')
linestyle -h h1 -4K 
linestyle -h h default
subplot(2,2,3),h=plot(r,zeta);vlines(10,'k:'),hlines(0,'k:'),ylim([-1 5.1]/1e5), xlim([0 50])
h1=h(1);h=h(2:end);title('Vorticity'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,4),h=plot(f./dr,Psi);xlog,ylog,ylim(10.^[-6 1]),vlines(1./10/2/pi,'k:'),xlim([1e-3 maxmax(f./dr)])
h1=h(1);h=h(2:end);title('Spectrum'),xlabel('Radian Wavenumber')
linestyle -h h1 -4K
linestyle -h h default
%letterlabels(2)

cd_figures
fontsize 14 12 12 12
orient landscape
print -depsc materneddies_core.eps


alpha=logspace(log10(3/4),log10(25),100);
figure,
plot(alpha,frac(maternfun(alpha-1,maternrad(alpha)),maternfun(alpha-1,0)))
linestyle -2k
hlines(1./sqrt(exp(1)),'k:')
title('Rossby Number Ratio'),xlabel('Shape Parameter \alpha'),xlim([0 25])
fontsize 14 12 12 12
orient portrait
%print -depsc rossbyratio.eps

%A la McWilliams
figure
h=plot((r./10).^2,abs(v));vlines(10,'k:'),hlines(0,'k:'),xlim([0 25]),ylog
h1=h(1);h=h(2:end);title('Velocity'),xlabel('Distance Squared'),ylim([.1 15])
orient portrait



function[]=morseddy_figures
dr=0.1;
r=[0:dr:500]';
alpha=[1.5:1/4:4];
%alpha=logspace(log10(1.25),log10(5),10);
[v,psi,zeta]=simpleddy(r,10,15,'morse',alpha);
[vo,psio,zetao]=simpleddy(r,10,15,'gaussian');


psi=[psio psi];
v=[vo v];
zeta=[zetao zeta];

psitwosides=[flipud(psi);psi(2:end,:)];
Psi=abs(fft(psitwosides));
%[mu,sigma,skew,kurt]=pdfprops(1:length(Psi),psitwosides(:,[2:end 1]));
%[muf,sigmaf,skewf,kurtf]=pdfprops(1:length(Psi),fftshift(Psi(:,[2:end 1])));

f=fourier(pi,length(Psi));
Psi=Psi(1:length(r),:);

figure
subplot(2,2,1),h=plot(r,psi);vlines(10,'k:'), xlim([0 50]),
h1=h(1);h=h(2:end);title('Streamfunction'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,2),h=plot(r,v);vlines(10,'k:'),ylim([0 15.5]), xlim([0 50])
h1=h(1);h=h(2:end);title('Velocity'),xlabel('Distance')
linestyle -h h1 -4K 
linestyle -h h default
subplot(2,2,3),h=plot(r,zeta);vlines(10,'k:'),hlines(0,'k:'),ylim([-3 10]/1e5), xlim([0 50])
h1=h(1);h=h(2:end);title('Vorticity'),xlabel('Distance')
linestyle -h h1 -4K
linestyle -h h default
subplot(2,2,4),h=plot(f./dr,Psi);xlog,ylog,ylim(10.^[-6 1]),vlines(1./10/2/pi,'k:'),xlim([1e-3 maxmax(f./dr)])
h1=h(1);h=h(2:end);title('Spectrum'),xlabel('Radian Wavenumber')
linestyle -h h1 -4K
linestyle -h h default
%letterlabels(2)

cd_figures
fontsize 14 12 12 12
orient landscape
%print -depsc morseddies.eps


% alpha=logspace(log10(3/4),log10(25),100);
% figure,
% plot(alpha,frac(maternfun(alpha-1,maternrad(alpha)),maternfun(alpha-1,0)))
% linestyle -2k
% hlines(1./sqrt(exp(1)),'k:')
% title('Rossby Number Ratio'),xlabel('Shape Parameter \gamma'),xlim([0 25])
% fontsize 14 12 12 12
% orient portrait
% print -depsc morserossbyratio.eps



function[]=eddyslice_figures

x=[-100:.1:100]';
y=[0:2:16];
clear cx
for i=1:length(y)
    cx(:,i)=x+sqrt(-1)*y(i);
end

Ro=10;
Vo=15;


[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'gaussian');
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',1/2);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',1);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo,'matern',2);
%[cv,psi,zeta]=simpleddy(cx,Ro,Vo);

figure,
subplot(2,2,1),plot(cx),hold on,axis([-1 1 -1 1]*Ro*2),axis equal
h=plot(Ro*rot(0:.01:2*pi+.1));
linestyle -h h 2G
subplot(2,2,2),hodograph(cv)
subplot(2,2,3),plot(x,real(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Along-Stream Currents')
subplot(2,2,4),plot(x,imag(cv)),vlines([-Ro Ro],'k:');hlines(0,'k:');
title('Cross-Stream Currents')

% 
% x=[-20:.1:20]';
% [xg,yg]=meshgrid(x,x);
% [cv,psi,zeta]=simpleddy(xg+sqrt(-1)*yg,Ro,Vo,'matern',1);




