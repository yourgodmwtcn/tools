% roms_obs_example.m
%
% Neil Banas, oct 2010
%
% this demonstrates the basic syntax for identifying a ROMS run,
% extracting some fields using roms_extract.m, and extracting
% some corresponding observations using obs_extract.m.



% what and where and when ------------------------------------------------------

RUN_DIR = '/pmraid1/daves/runs/salish_2006_3/OUT/';
OBS_DIR = '/skua1/neil/mossea_obs_data/';
SECTION_DIR = '/skua1/neil/mossea_geo_data/sections/';

% locate the _his files (run.his) and the model grid info (run.grid)
run = roms_createRunDef('salish_2006_3', RUN_DIR);

% identify the low-passed files too, in an adjacent directory
run.lp = roms_createSeriesDef([RUN_DIR '../OUT_lp/'], 'ocean_his_lp_');

% timebase: model save times from Jul 8-18, 2006
T = run.his.nctime(run.his.nctime >= datenum('7/8/2006') & ...
                   run.his.nctime <= datenum('7/18/2006'));
T = T(1:6:end); % to speed things up

% a section line and a station along it
load([SECTION_DIR 'JdF_PS_section.mat']); % creates 'Section'
whichLabel = strmatch('Main Basin',Section.label);
ind = Section.label_ind(whichLabel);
Stn.x = Section.x(ind);
Stn.y = Section.y(ind);



% ------------------------------------------------------------------------------
% plot some coordinated extractions of model data. I've omitted all the
% graphical niceities.
figure;
contourLevels = [0:4:28 28.5:0.5:33];

subplot 221
% extract surface salinity at one time point and plot it
[S, coords] = roms_extract(run.his, 'salt', T(1), 'surface');
contourf(coords.xm, coords.ym, S, contourLevels, 'edgecolor', 'none');
title('surface salinity on Jul 8');
hold on;
plot(Section.x, Section.y, 'k');
plot(Stn.x, Stn.y, 'k*');

subplot 222
% extract a salinity section at one time
[S, coords] = roms_extract(run.his, 'salt', T(1), 'profile', Section.y, Section.x);
% make an along-section distance variable out of the coordinate variables coords.xm, coords.ym
coords.dist = trackDist(coords.xm(:), coords.ym(:), Section);
coords.dist = reshape(coords.dist, size(S));
% plot it
contourf(coords.dist, coords.zm, S, contourLevels, 'edgecolor', 'none');
title('salinity from JdF - Main Basin, Jul 8');

subplot 223
% depth vs time at the station
[S, coords] = roms_extract(run.his, 'salt', T, 'profile', Stn.y, Stn.x);
% to do this for low-passed salinity, it's just 
% [S, coords] = roms_extract(run.lp, 'salt', T, 'profile', Stn.y, Stn.x);
contourf(coords.tm, coords.zm, S);
datetick('x');
title('salinity at Main Basin station, Jul 8-18');


% ------------------------------------------------------------------------------
% find all CTDs within 5 km of the section line during this period
obs = obs_extract(OBS_DIR, 'salinity', T([1 end]), 'section', Section.x, Section.y, 5);

% in this example the following lines are unnecessary,
% but you could clean up bad points with filters like
obs = obs_omit(obs, obs.salinity < 0 | obs.salinity > 40);
% or select subsets with filters like
obs_part1 = obs_omit(obs, obs.t > mean(T));
obs_part2 = obs_omit(obs, obs.t <= mean(T));

figure;
subplot 121
% scatter all CTD salinity vs depth 
plot(obs.salinity,-obs.z,'k.');
title('salinity vs depth, CTDs JdF-Main Basin, Jul 8-18');

subplot 122;
% to look at the data cast by cast...
casts = obs_separateCasts(obs);
% now plot each cast against the model at the same x,y,t
for i=1:length(casts)
	disp([num2str(i) ' out of ' num2str(length(casts))]);
	modelS = roms_extract(run.his, 'salt', casts(i).t(1), 'point', -casts(i).z, casts(i).y, casts(i).x);
	plot(casts(i).salinity, modelS, 'k');
	hold on;
end
plot([23 33],[23 33],'k:');
xlabel('CTD');
ylabel('model');
title('salinity from JdF - Main Basin, Jul 8-18 (CTD vs model)');
