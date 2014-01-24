function []=filter_test(M)

% filter_test.m 8/19/2013 Parker MacCready
%
% plots the results of a mooring extraction and then applies some
% filters to

td = M.td;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

figure; set(gcf,'position',[20 20 1400 900]); Z_fig;

eta = M.zeta; 

subplot(411)
plot(td0,eta)
ylabel('Sfc. Height (m)')
grid on
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')
xlim([td0(1) td0(end)]);
aa = axis;

subplot(412)
plot(td0,Z_jfilt(eta',40))
title('40 hour Hanning Window');
grid on
xlim([td0(1) td0(end)]);
aa = axis;

subplot(413)
plot(td0,Z_dasfilt(eta,'godin'))
title('Godin 24-24-25 Filter');
grid on
axis(aa);

subplot(414)
plot(td0,Z_dasfilt(eta,'ttide'))
title('T Tide Filter');
grid on
axis(aa);

xlabel('Yearday')


