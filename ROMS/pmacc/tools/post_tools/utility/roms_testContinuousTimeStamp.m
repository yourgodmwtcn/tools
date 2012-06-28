function [t, dt] = roms_testContinuousTimeStamp(dirname, ncbasename)

%this function tests that the files within dirname (which should be the
%ROMS history output directory, or could be a lp directory) have a
%continuos uniform time interval
%
%dirname is the directory name such as /pmraid2/sarahgid/runs/runname/OUT/
%
%ncbasename is the basename of interest, usually ocean_his_ but could also
%be ocean_his_lp_ or some other name of choice
%
%SNG April 2011

addpath(genpath('/home/sarahgid/mexcdf/'))

if nargin < 2
    ncbasename = 'ocean_his_';
end

filelist = dir([dirname ncbasename '*.nc']);
ncnames = {filelist(:).name};
% extract the file numbers
ncnums = strrep(ncnames, ncbasename, '');
ncnums = strrep(ncnums, '.nc', '');
ncn = nan.*ones(1,length(ncnums));
for i=1:length(ncnums)
	n = str2double(ncnums{i});
	if length(n) ~= 1, n = nan; end
	ncn(i) = n;
end
ncnames = ncnames(~isnan(ncn));

t = nan.*ones(1,length(ncnums));
for i = 1:length(ncnums)
    filein = ncnames{i};
    t(i) = nc_varget([dirname filein],'ocean_time');
end

dt = diff(t);
figure;
plot(dt);
if isempty(find(dt~=dt(1), 1))
    disp([dirname ' contains history files with a uniform time interval of ' num2str(dt(1)./(60*60)) ' hours'])
else
    disp(['WARNING! ' dirname ' does NOT contain a uniform time interval!'])
end