function [h]=timeplt(jd,u,istack,ylims);
% TIMEPLT  time series stack plots with Gregorian Time labels on x axes
%   
%      type timeplt('demo') for a demonstration or
%
%      USAGE:  
%          [h]=timeplt(jd,u,[istack],[ylims]);
%      where
%         jd = Julian Day time vector  (e.g produced by JULIAN.M)
%         u = column vector or matrix of column vectors containing time
%             series data.  If the column is complex, it will be plotted
%             as a stick plot.  
%         istack = vector of indices indicating which panel you want
%             to plot the time series data.  istack=[1 2] would make
%             two panels one on top of the other and plot the first 
%             column of u in the lower panel and the second column of
%             u in the upper panel.  If any column in u is complex,  
%             istack must be specified.  If istack is not specified, all the 
%             columns will be plotted in the first panel.
%         ylims  = [npanels x 2] matrix containing the ylimits of 
%              the panel plots.  If you are plotting two panels and 
%              you want the limits of both plots to be from -10 to 15,
%              then set ylims=[-10 15; -10 15].  Autoscales if ylims
%              is not set
%
%       outputs:
%            h = handles for stack plots (axes)
%
%
% requires GREGAXD, GREGAXM, GREGAXH, GREGAXY, GREGORIAN, JULIAN, and
%          STACKLBL
%
% Rich Signell rsignell@usgs.gov
%
% 8-29-94 fixed bug in stick orientation 
% 8-30-94 fixed problem with upper panel too close to top of figure
% 11-1-94 eliminated ylabels from the argument list.  Just use the
%         axes handles returned by timeplt to change axes and then
%         use the standard "ylabel" command.
% 12-4-96 replaced vector stick plotting method from lots of 2 point
%         objects to 1 object with nans.  Much faster plotting. 
%

%-------------------------------------------------------------------
%  Set the cutoff for different types of Gregorian axis types
%  You can adjust these to suit your preferences.  For example, if your
%  plot got labeled with days, but you want hours, increase "daycut" 
%  (until it's larger than the "fac" for your plot).
  
  yearcut=250;
  moncut=20;
  daycut=.2;
  mincut=0.05;
%
if strcmp(jd,'demo')
  help timeplt
  start=[1990 11 1 0 0 0];    %Gregorian start [yyyy mm dd hh mi sc]
  stop=[1991 2 1 0 0 0];
  jd=julian(start):julian(stop); 
  u=sin(.1*jd(:)).^2-.5;
  v=cos(.1*jd(:));
  w=u+i*v;
%w is vector, so must have it's own axes
  h=timeplt(jd,[u v abs(w) w],[1 1 2 3]); 
  title('Demo of Timeplt')
  stacklbl(h(1),'East + North velocity','m/s');
  stacklbl(h(2),'Speed','m/s');
  stacklbl(h(3),'Velocity Sticks','m/s');
  return
end


%
% Don't clip that data!
set(0,'DefaultLineclipping','off')

[m,n]=size(u);
if(m==1),u=u(:);m=n;n=1;end; %if row vector, make column vec
clf
if(~exist('istack')),
  istack=ones(n,1);
end
njd=length(jd);
jd0=jd(1);
jd1=jd(njd);
ndays=jd1-jd0;
%
% determine how many stack plots to make
%
nstack=max(istack);
% determine height of each stack plot
hstack=.70/nstack;
xmin=.13;
ymin=.15;
width=0.77;
h(1)=axes('position',[xmin ymin width hstack]);
set(gca,'units','pixels');...
pos=get(gca,'pos');
xlen=pos(3);
ylen=pos(4);...
font=get(gca,'fontsize');...
label_width=5*font;
nlabel=floor(xlen/label_width);
set(gca,'units','normalized');
fac=ndays/nlabel;  %number of days per label
% adjust xfactor for subsequent stretching by Gregorian Date
  if( fac > yearcut),
    xlim=[jd0-180 jd1+180];
  elseif(yearcut > fac & fac > moncut),
    xlim=[jd0-15 jd1+15];
  elseif(moncut > fac & fac > daycut),
    xlim=[floor(jd0)-.5 ceil(jd1)+.5];
  elseif(daycut > fac & fac > mincut);
    xlim=[jd0-1/48 jd1+1/48];
  elseif(mincut > fac);
	
	%
	% leave 2.5% off of each side.
	time_offset = (jd1-jd0)*0.025;
    xlim=[jd0-time_offset jd1+time_offset];
  end
up=u(:,find(istack==1));
if(isreal(up)),
  h2=plot(jd,up);... 
  if(exist('ylims')),
    set(gca,'ylim',ylims(1,:));...
  end
else
  if(exist('ylims')),
    ylim=ylims(1,:);
  else
    y0=min(gmin(imag(up(:))),0);
    y1=max(gmax(imag(up(:))),0);
    ylim=[y0 y1];
  end
  set(gca,'units','pixels');
  ppos=get(gca,'position');
  set(gca,'units','norm');
  d=diff(xlim);
  uscale=(diff(xlim)/diff(ylim))*(ppos(4)/ppos(3));
  vp=imag(up);
  up=uscale.*real(up);
  x=jd;
  xp=x;
  yp=zeros(size(xp));
  xplot=ones(length(xp),2);
  yplot=xplot;
  xplot(:,1)=x(:);
  xplot(:,2)=xp(:)+up(:);
  xplot(:,3)=x(:);
  yplot(:,1)=yp(:);
  yplot(:,2)=yp(:)+vp(:);
  yplot(:,3)=yp(:)*nan;
  xplot=xplot';
  yplot=yplot';
  if(~isempty(find(finite(up(:))))),
   plot([jd0 jd1],[0 0],'r-',xplot(:),yplot(:),'r-');...
  set(gca,'ylim',ylim);
  end
end
set(gca,'xlim',xlim);
if( fac > yearcut),
  gregaxy(jd,floor(fac/yearcut));
elseif(yearcut >= fac & fac>moncut),
  gregaxm(jd,floor(fac/moncut));
elseif(moncut >= fac & fac>daycut),
  gregaxd(jd,ceil(fac));
elseif(daycut >= fac & fac> mincut);
  gregaxh(jd,max(1,floor(fac*48)));
elseif(mincut >= fac);
  gregaxmi(jd,max(1,floor(fac*120*24)));
end
if(nstack==1),return,end
pos_norm=get(gca,'position');
vsep=1.3*font*pos_norm(4)/pos(4);
for iplot=2:nstack,
  h(iplot)=axes('position',[xmin ymin+(hstack+vsep)*(iplot-1) width hstack]);
  up=u(:,find(istack==iplot));
if(isreal(up)),
  plot(jd,up);
  if(exist('ylims')),
    set(gca,'ylim',ylims(iplot,:))
  end
else
  if(exist('ylims')),
    ylim=ylims(iplot,:);
  else
    y0=min(0,gmin(imag(up(:))));
    y1=max(0,gmax(imag(up(:))));
    ylim=[y0 y1];
  end
  set(gca,'units','pixels');
  ppos=get(gca,'position');
  set(gca,'units','norm');
  d=diff(xlim);
  uscale=(diff(xlim)/diff(ylim))*(ppos(4)/ppos(3));
  vp=imag(up);
  up=uscale.*real(up);
  x=jd;
  xp=x;
  yp=zeros(size(xp));
  xplot=ones(length(xp),2);
  yplot=xplot;
  xplot(:,1)=x(:);
  xplot(:,2)=xp(:)+up(:);
  xplot(:,3)=x(:);
  yplot(:,1)=yp(:);
  yplot(:,2)=yp(:)+vp(:);
  yplot(:,3)=yp(:)*nan;
  xplot=xplot';
  yplot=yplot';
  if(~isempty(find(finite(up(:))))),
  plot([jd0 jd1],[0 0],'r-',xplot(:),yplot(:),'r-');...
  set(gca,'ylim',ylim);
  end
end
  set(gca,'xlim',xlim,'xtick',get(h(1),'xtick'),'xticklabels',[]);
end
set(0,'DefaultLineclipping','on')
