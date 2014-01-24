%data = obs_extract(files, vars, timeRange, varargin)

%% files

files='C:\Documents and Settings\Corinne\MatlabWork\PSVS\'

%% vars
vars={'salinity','temperature'}
%% timerange
%% polygon
vars={'salinity','temperature','oxygen','ob1'};
load('C:\Documents and Settings\Corinne\MatlabWork\PSVS\HoodCanalModeling\Thalweg_5.mat');
 ring=make_range_ring(lon(1), lat(1),3);%lon,;lat,,range
files ='C:\Documents and Settings\Corinne\MatlabWork\PSVS\netcdf_files';
data = obs_extract(files, vars, [datenum(2004,1,1) datenum(2008,1,1)], 'polygon',ring(1),ring(2));

%% section

load('C:\Documents and Settings\Corinne\MatlabWork\PSVS\HoodCanalModeling\Thalweg_5.mat');
files ='C:\Documents and Settings\Corinne\MatlabWork\PSVS\';
data = obs_extract(files, {'salinity','temperature','oxygen','nitrate_bottle'}, [datenum(2006,1,1) datenum(2007,1,1)], 'section',lon,lat,3);
 plot(data.oxygen(data.z>=60),data.temperature(data.z>=60),'.')
plot(data.salinity(data.z>=60),data.nitrate_bottle(data.z>=60),'.')

%% one point