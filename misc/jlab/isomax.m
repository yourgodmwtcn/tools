function[varargout]=isomax(varargin)
% ISOMAX  Locates isolated maximum of wavelet spectrum.
%  
%   [IX,JX]=ISOMAX(W,CHI) returns the locations (II,JJ) of all isolated
%   maxima of the wavelet transform W having magnitudes exceeding CHI, i.e. 
%   ABS(W)>CHI.  IX and JX are column vectors with IX(m) and JX(m) 
%   specifying the row and column, respectively, of the mth isolated maxima.
%   The maxima are sorted in the order of decreasing strength. 
%
%   [IX,JX]=ISOMAX(Y,CHI,DT,DF,FS) only finds those maxima which are the 
%   are the largest magnitude points within a time frequency box of time 
%   half-width DT and frequency half-width DF. Here FS is an array
%   specifying the frequencies corresponding to the columns of Y.  DT and 
%   DF may each be scalars or arrays of length LENGTH(FS).
%
%   [IX,JX]=ISOMAX(Y,CHI,DT,DF,FS,SMIN) only finds those maxmima occuring 
%   at scale numbers (e.g. column numbers) greater than SMIN.
%
%   [IX,JX,WX]=ISOMAX(...) also returns the corresponding value of the 
%   wavelet transform W at each isolated maximum.
%
%   Usage: [ix,jx,wx]=isomax(y,chi);
%          [ix,jx,wx]=isomax(y,chi,dt,df,fs);
%          [ix,jx,wx]=isomax(y,chi,dt,df,fs,smin);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2012 J.M. Lilly --- type 'help jlab_license' for details        

x=varargin{1};
a=varargin{2};
dt=[];
df=[];
wf=[];
smin=[];
if nargin>=3
  dt=varargin{3};
end
if nargin>=4
  df=varargin{4};
end
if nargin>=5
  wf=varargin{5};
end
if nargin>=6
  smin=varargin{6};
end



for i=1:size(x,3)
   [ii{i},jj{i},xm{i}]=isomax1(x(:,:,i),a,dt,df,wf,smin);
   kk{i}=i+zeros(size(ii{i}));
end
vcellcat(ii,jj,kk,xm);
varargout{1}=ii;
varargout{2}=jj;
if size(x,3)==1
  varargout{3}=xm;
else
  varargout{3}=kk;
  varargout{4}=xm;
end

function[ii,jj,xm]=isomax1(x,a,dt,df,wf,smin)
x=abs(x);

%Check magnitude of X versus its neighbors  
b1=x>vshift(x,1,1);
b2=x>vshift(x,-1,1);
b3=x>vshift(x,1,2);
b4=x>vshift(x,-1,2);
bool=b1.*b2.*b3.*b4;

%Exclude points on edges
bool(1,:)=0;
bool(end,:)=0;
bool(:,1)=0;
bool(:,end)=0;

if ~isempty(smin)
  bool(:,1:smin)=0;
end

%Find value of X at isolated maxima
index=find(bool);
xm=x(index);

%Remove points below cutoff
index2=find(xm>a);
if ~isempty(index2)
  vindex(index,xm,index2,1);
  
  %Sort remainer
  [temp,index3]=sort(-xm);
  xm=xm(index3);
  index=index(index3);
  
  [ii1,jj1]=ind2sub(size(x),index);
  bool=1+0*ii1;
  
  if isempty(dt)
      ii=ii1;
      jj=jj1;
  else
      %Search within boxes 
      wf=wf(:);
   
      if length(dt)==1
            dt=dt+0*fs;
      end
      if length(df)==1
            df=df+0*fs;
      end
      
        
      for i=1:length(xm)
        imin=ii1(i)-dt(jj1(i));
        imax=ii1(i)+dt(jj1(i));
        
        fmin=wf(jj1(i))-df(jj1(i));
        fmax=wf(jj1(i))+df(jj1(i));
        
        %Locate those within same box but smaller amplitude
        %ii1,imax,imin, xm
        
        b1=(ii1<imax & ii1>imin) & xm<xm(i);
        b2=(wf(jj1)<fmax & wf(jj1)>fmin) & xm<xm(i);
        index=find(b1 & b2);
        %b2count(i)=length(find(b2));
        %b1count(i)=length(find(b1));
        
        %These are excluded from census
        if ~isempty(index)
          bool(index)=0;
        end
     end
     %Extract those which have not been excluded
     index=find(bool);
     if ~isempty(index)
        [ii,jj,xm]=vindex(ii1,jj1,xm,index,1);
     end
  end
end

function[]=isomax_test;

load bravo94
use bravo.rcm


%Decide on frequencies
fs=morsespace(1,2,2*pi/50,2*pi/1000,8);

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cv(:,2)),imag(cv(:,2)),{1,1,2,fs,'bandpass'},'mirror');
w=sqrt(abs(wx).^2+abs(wy).^2);
figure,
pcolor(w'),shading interp,flipy,hold on
[ix,jx,wx]=isomax(w,10,2./fs,fs/2,fs);
plot(ix,jx,'m*')


