% coast_plot.m 5/5/2011 Parker MacCready
%
% plot the results of coast_mat_maker.m

% calls code from pandora/Z_functions

clear

c1 = load('pnw_coast_detailed.mat');
c2 = load('pnw_coast_regional.mat');
c3 = load('pnw_coast_combined.mat');
c4 = load('pnw_rivers.mat');

clf
Z_fig(10)

subplot(131)
plot(c1.lon,c1.lat,'-r');
Z_dar;
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title('\color{red}detailed')

subplot(132)
plot(c2.lon,c2.lat,'-k');
Z_dar;
hold on
plot(c4.lon,c4.lat,'-g');
xlabel('Longitude (deg)')
title('regional + \color{green}rivers')

subplot(133)
plot(c3.lon,c3.lat,'-b');
Z_dar;
xlabel('Longitude (deg)')
title('\color{blue}combined')

set(gcf,'position',[100 100 1000 600]);
set(gcf,'PaperPositionMode','auto');
print('-djpeg90','coast_plot.jpg')
