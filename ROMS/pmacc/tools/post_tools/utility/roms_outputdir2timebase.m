function [nctime, ncn] = roms_outputdir2timebase(dirname,ncbasename)

% [nctime, ncn] = roms_outputdir2timebase(dirname,ncbasename);
%
% given a directory containing a bunch of single-frame ROMS output netcdf files,
% returns a pair of variables that associate filenumbers (ncn) with model times (nctime).
%
% assumes a uniform timebase: only looks inside the first and last file to construct nctime.
%
% example:
% [nctime, ncn] = outputdir2timebase('~/cr5_2004_bio9/OUT/', 'ocean_his_');
%
% neil banas feb 2009

ncn = [];
nctime = [];
if exist(dirname) ~= 7, return; end

% locate files in the directory _dirname_ that look like 'ncbasename*.nc'
filelist = dir([dirname ncbasename '*.nc']);
ncnames = {filelist(:).name};
% extract the file numbers
ncnums = strrep(ncnames, ncbasename, '');
ncnums = strrep(ncnums, '.nc', '');
for i=1:length(ncnums)
	n = str2num(ncnums{i});
	if length(n) ~= 1, n = nan; end
	ncn(i) = n;
end
ncnames = {ncnames{~isnan(ncn)}};
ncn = ncn(~isnan(ncn));

if isempty(ncn), return; end
% look in the first and the last files to get the first and last model times
nctime = repmat(nan, size(ncn));
nctime(1) = nc_varget([dirname ncnames{1}],'ocean_time');
nctime(end) = nc_varget([dirname ncnames{end}],'ocean_time');
units = nc_attget([dirname ncnames{1}],'ocean_time','units');
timeref = strrep(units, 'seconds since ', '');

% interpolate times
nctime(2:end-1) = interp1(ncn([1 end]), nctime([1 end]), ncn(2:end-1));

% convert nctime to matlab date format
nctime = nctime./86400 + datenum(timeref,'yyyy-mm-dd HH:MM:SS');

disp(['found ' num2str(length(ncn)) ' files from ' ncnames{1} ' to ' ncnames{end} '; ']);
disp(['times from ' datestr(nctime(1)) ' to ' datestr(nctime(end))]);