function[v]=ellvel(varargin)
% ELLVEL  Average and instantaneous ellipse velocities. 
%
%   V=ELLVEL(KAPPA,LAMBDA,THETA,PHI,STR) where KAPPA and LAMBDA are the 
%   amplitude and linearity of a time-varying ellise, THETA is its time-
%   varying orientation, and PHI is its time-varying phase, returns various
%   measures of the ellispe 'velocity'.
%
%   STR determines the velocity quantity to be output.
%
%       STR       Symbol     Description
%       ----------------------------------------------------------------
%      'geo'      VM         Geometric mean velocity    
%      'ave'      VBAR       Period-averaged speed          
%      'kin'      VEKE       Kinetric energy velocity        
%      'cir'      VGAMMA     Circulation velocity        
%      'azi'      VAZ        Instantaneous azimuthal velocity        
%      'ins'      VI         Instantaneous speed 
%
%   All are signed quantities reflecting the direction of motion.  That is, 
%   a velocity is defined to be positive when the ellipse is orbited in the 
%   mathematically positive (counterclockwise) sense, and negative when the
%   ellipse is orbited in the mathematically negative sense.
%
%   See Lilly and Gascard (2006) for details on all quantities except for 
%   the circulation velocity VGAMMA and VEKE, which are discussed in XX.
%
%   VM=ELLVEL(KAPPA,LAMBDA,THETA,PHI) with STR omitted returns VM, the
%   geometric mean velocity.  This is a basic of ellipse velocity that is 
%   analagous to the geometric mean radius as a measure of ellipse size.
%
%   ELLVEL(DT,...,) optionally uses DT as the data sample rate, with a 
%   default value of DT=1.  DT is a scalar.
%
%   ELLVEL(DT,...,FACT,STR) optionally converts the physical units of 
%   velocity through a multiplication by FACT, with a default value of 
%   FACT=1.  For example, FACT=1e5 converts kilometers into centimeters.
%   ____________________________________________________________________
%   
%   Cell array input/output
%
%   If ELLVEL is given cell array input, it returns cell array output.
%
%   Thus KAPPA, LAMBDA, THETA, and PHI may each be cell arrays of the 
%   same size, where each element in the cell array is a numerical array.
%
%   VM and other velocity measures will be also cell arrays of this size.
%
%   In this case ELLVEL(DT,...) also works with DT a scalar or an 
%   array whose length is the number of elements in the cell arrays.
%   ____________________________________________________________________
%
%   See also ELLRAD, ELLPARAMS, ELLDIFF.
%
%   Usage:  v=ellvel(kappa,lambda,theta,phi);
%           v=ellvel(dt,kappa,lambda,theta,phi,1e5);
%           vgamma=ellvel(kappa,lambda,theta,phi,'circulation');
%           vgamma=ellvel(dt,kappa,lambda,theta,phi,fact,'circulation');
%
%   'ellvel --f' generates a sample figure.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2012 J.M. Lilly --- type 'help jlab_license' for details    

%       VPHIBAR    Period-averaged azimuthal velocity 

if strcmp(varargin{1}, '--f')
  ellvel_fig,return
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else 
    dt=1;
end

if isstr(varargin{end})
    str=varargin{end}(1:3);
    varargin=varargin(1:end-1);
else
    str='geo';
end

if length(varargin{end})==1
    fact=varargin{end};
    varargin=varargin(1:end-1);
else 
    fact=1;
end

kappa=varargin{1};
lambda=varargin{2};
theta=varargin{3};
phi=varargin{4};

if ~iscell(kappa)
    v=ellvel_one(dt,kappa,lambda,theta,phi,fact,str);
else
    for i=1:length(kappa)
%         if isscalar(dt)
%             dti=dt;
%         else
%             dti=dt(i);
%         end
%         if isscalar(fact)
%             facti=fact;
%         else 
%             facti=fact{i};
%         end
        v{i}=ellvel_one(dt,kappa{i},lambda{i},theta{i},phi{i},fact,str);
    end
end




function[v]=ellvel_one(dt,kappa,lambda,theta,phi,fact,str)

strend='endpoint';
vswap(lambda,0,1e-10);
vswap(lambda,1,1-1e-10);

if aresame(str,'geo')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2);
elseif aresame(str,'ave')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2,'ave');
elseif aresame(str,'cir')
    omphi=frac(1,dt).*vdiff(unwrap(angle(rot(phi))),1,strend);
    omtheta=frac(1,dt).*vdiff(unwrap(angle(rot(theta))),1,strend);
    om=omphi+sign(lambda).*sqrt(1-lambda.^2).*omtheta;
    v=om.*kappa.*sqrt(1-lambda.^2);
elseif aresame(str,'kin')
    v=elldiff(dt,kappa,lambda,theta,phi);
elseif aresame(str,'azi')
    z=ellsig(kappa,lambda,theta,phi);
    Omega=frac(1,dt).*vdiff(unwrap(angle(z)),1,str);
    v=Omega.*ellrad(kappa,lambda,phi,'instantaneous');
elseif aresame(str,'ins')
    [k2,l2,theta2,phi2]=elldiff(dt,kappa,lambda,theta,phi);
    v=ellrad(k2,l2,phi2,'instantaneous');
end

v=abs(v).*sign(lambda).*fact;    


function[]=ellvel_test

load ebasnfloats
use ebasnfloats


len=cellength(lat);
index=find(len>200);
lato=30;

id=id(index);num=num(index);lat=lat(index);lon=lon(index);
p=p(index);t=t(index);

index=24;
id=id(index);num=num{index};lat=lat{index};lon=lon{index};
p=p{index};t=t{index};
dt=num(2)-num(1);

ga=3;be=3;
fs=morsespace(ga,be,{0.05,2*pi/3},2*pi/100,8);


%Compute wavelet transforms using generalized Morse wavelets

cx=fillbad(latlon2xy(lat,lon,30,-25));
cv=latlon2uv(num,lat,lon);

wx=wavetrans(real(cx),{1,ga,be,fs,'bandpass'},'mirror');
wy=wavetrans(imag(cx),{1,ga,be,fs,'bandpass'},'mirror');

[ir,jr,xr,fr]=ridgewalk(dt,wx,wy,fs,{2*morseprops(ga,be),0,'amp'});  
[xhat,fhat]=ridgemap([length(cx) 2],xr,fr,ir);
fbar=instfreq(dt,xhat,2);
[kap,lam,the,phi]=ellparams(xhat(:,1),xhat(:,2));  

R=ellrad(kap,lam,phi);
[vm,vgamma]=ellvel(dt,kap,lam,the,phi);
%This test is not done!! 

function[]=ellvel_fig
lambda=(0:.001:1)';
kappa=1+0*lambda;
phi=(1:length(lambda))'/10;
theta=0*lambda;

fact=1e5;dt=1/24;
vm=ellvel(dt,kappa,lambda,theta,phi,fact,'geometric');
vgamma=ellvel(dt,kappa,lambda,theta,phi,fact,'circulation');
veke=ellvel(dt,kappa,lambda,theta,phi,fact,'kineticenergy');
vbar=ellvel(dt,kappa,lambda,theta,phi,fact,'average');

figure
plot(lambda,[vgamma,vm,vbar,veke]./maxmax(vm));
linestyle 2k k-- k 2G
e=[.2486 .967];
axis([0.01 1 0 1.05]),axis square 
vlines(e.^2./(2-e.^2),'k:')
%axis([0 1 0.5 1]),axis square 
legend('V_\Gamma','V_M','V_{Bar}','V_{EKE}')
title('Velocity measures for constant \kappa and frequency')
xlabel('Ellipse linearity \lambda')
ylabel('Mean velocity measures')
set(gcf,'paperposition', [2 2 3.5 3.5])
xtick(.1),ytick(.1),fixlabels(-1)
%text(0.75,0.93,'V_{Bar}')
%text(0.55,0.83,'V_M') 
fontsize 14 14 14 14
%fontsize jpofigure
%cd_figures
%print -depsc ellipsemeans.eps


