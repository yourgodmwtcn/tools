function []=basic(M)

% basic.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction

td = M.td;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

figure; set(gcf,'position',[20 20 1400 900]); Z_fig;

subplot(411)
plot(td0,M.zeta)
ylabel('Sfc. Height (m)')
grid on
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')
xlim([td0(1) td0(end)]);

subplot(412)
plot(td0,M.sustr,'-b',td0,M.svstr,'-r')
legend('sustr','svstr',0)
ylabel('Wind Stress (Pa)')
grid on
xlim([td0(1) td0(end)]);


subplot(413)
plot(td0,M.shflux)
ylabel('shflux (W)')
grid on
xlim([td0(1) td0(end)]);


subplot(414)
plot(td0,M.u(end,:),'-b',td0,M.v(end,:),'-r')
legend('u','v',0)
xlabel('Date')
ylabel('Sfc. Vel. (m s^{-1})')
grid on
aa = axis;
xlim([td0(1) td0(end)]);

