function []=hflux(M)

% hflux.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction, focused on surface heat flux

td = M.td;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

td00 = td0(1); td11 = td0(end);
mask = find(td0>=td00 & td0<=td11);
td0 = td0(mask);

ymin = -500; ymax = 500;

figure; set(gcf,'position',[20 20 1400 900]); Z_fig(16);

subplot(411)
plot(td0,M.shflux(mask))
ylabel('shflux')
grid on
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')
axis([td00 td11 ymin ymax]);
[xt,yt] = Z_lab('ul');
text(xt,yt,['Mean = ',num2str(round(mean(M.shflux(mask)))), ...
    ' (W m^{-2})'],'fontweight','bold');

subplot(412)
plot(td0,M.latent(mask),'-b',td0,M.sensible(mask),'-r')
legend('latent','sensible',0)
grid on
axis([td00 td11 ymin ymax]);
[xt,yt] = Z_lab('ul');
text(xt,yt,['Latent Mean = ',num2str(round(mean(M.latent(mask)))), ...
    ' (W m^{-2})'],'fontweight','bold');
[xt,yt] = Z_lab('ll');
text(xt,yt,['Sensible Mean = ',num2str(round(mean(M.sensible(mask)))), ...
    ' (W m^{-2})'],'fontweight','bold');

subplot(413)
plot(td0,M.swrad(mask))
ylabel('swrad')
grid on
axis([td00 td11 ymin ymax]);
[xt,yt] = Z_lab('ul');
text(xt,yt,['Mean = ',num2str(round(mean(M.swrad(mask)))), ...
    ' (W m^{-2})'],'fontweight','bold');


subplot(414)
plot(td0,M.lwrad(mask),'-r')
xlabel('Yearday')
ylabel('lwrad')
grid on
axis([td00 td11 ymin ymax]);
[xt,yt] = Z_lab('ul');
text(xt,yt,['Mean = ',num2str(round(mean(M.lwrad(mask)))),' (W m^{-2})'], ...
    'fontweight','bold');

