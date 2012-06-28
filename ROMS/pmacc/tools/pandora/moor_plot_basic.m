% moor_plot_basic.m 3/20/2012 Parker MacCready
%
% plots the results of a mooring extraction

clear;
[Tdir] = pan_start;
indir = [Tdir.pan_results,'moorings/'];
[fn,pth]=uigetfile([indir,'*.mat'], ...
    'Select mooring file...');
load([pth,fn]);
td = mod_moor.run.his.nctime;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

figure; set(gcf,'position',[20 20 1400 900]); Z_fig;

subplot(411)
plot(td0,mod_moor.zeta)
ylabel('Sfc. Height (m)')
grid on
name = strrep(fn,'_','\_'); name = strrep(name,'.mat','');
title(name,'fontweight','bold')
axis([td0(1) td0(end) -3 5]);

subplot(412)
plot(td0,mod_moor.sustr,'-b',td0,mod_moor.svstr,'-r')
legend('sustr','svstr',0)
ylabel('Wind Stress (Pa)')
grid on
axis([td0(1) td0(end) -.5 .5]);


subplot(413)
plot(td0,mod_moor.shflux)
ylabel('shflux (W)')
grid on
axis([td0(1) td0(end) -200 800]);


subplot(414)
plot(td0,mod_moor.u(end,:),'-b',td0,mod_moor.v(end,:),'-r')
legend('u','v',0)
xlabel('Date')
ylabel('Sfc. Vel. (m s^{-1})')
grid on
aa = axis;
axis([td0(1) td0(end) aa(3) aa(4)]);

