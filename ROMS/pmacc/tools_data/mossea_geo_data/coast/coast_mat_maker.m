% coast_mat_maker.m 5/2/2011 Parker MacCready
%
% this makes mat files out of a file created by the
% USGS Coastline Extractor

clear
close all

% Detailed (NOAA/NOS Medium Resolution Coastline) No Vancouver Island!
clear
a = load('pnw_coast_detailed.txt');
lon = a(:,1);
lat = a(:,2);
% add points that define the Columbia River mouth south jetty
lat = [lat; NaN; 46.23; 46.235];
lon = [lon; NaN; -124.02; -124.07];
save('pnw_coast_detailed.mat','lon','lat');

% Regional (World Vector Coastline)
clear
a = load('pnw_coast_regional.txt');
lon = a(:,1);
lat = a(:,2);
save('pnw_coast_regional.mat','lon','lat');

% Combined
clear
c1 = load('pnw_coast_detailed.mat');
c2 = load('pnw_coast_regional.mat');
lat = [c1.lat; NaN; c2.lat];
lon = [c1.lon; NaN; c2.lon];
save('pnw_coast_combined.mat','lon','lat');

% Rivers
clear
a = load('pnw_rivers.dat');
lon = a(:,1);
lat = a(:,2);
save('pnw_rivers.mat','lon','lat');

