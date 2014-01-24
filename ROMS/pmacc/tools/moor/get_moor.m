% get_moor.m  8/20/2013  Parker MacCready
%
% This is for extracting the results of ROMS simulations, focused on
% time series at a point, over the full water column, as would be
% collected by a mooring.
%
% NOTE: this calls the extraction code:
% [M] = Z_get_moor(run,mlon,mlat);
% which uses the native matlab netcdf calls, and so is relatively fast.

clear; close all; moor_start;

% &&&&&&&&&&& USER EDIT THIS &&&&&&&&&&&&&&
% First define the directory where the ROMS output files are.
% The "basename" is just a useful way of defining a run, and we use it in
% the code only for naming the "outfile"
basename = 'D2005';
indir = ['/Users/PM3/Documents/roms/output/',basename,'/'];
% The code allows you to extract from a bunch of mooring locations at once,
% and the cell array "mloc_list" is where you give names to each mooring.
% E.g. mloc_list = {'name1','name2'}; would allow two moorings.
mloc_list = {'RN'};
% The for-loop below automatically creates the vectors "mlon" and "mlat"
% which are contain the positions of the morings to be extracted.  These
% should be the same length as mloc_list.
% In this case we use the function "Z_mooringLocations" where we have
% stored information about moorings where we have observations to compare
% to, however there is no need to use this function.  All that is required
% is that you specify mlon and mlat vectors that define positions within
% the domain of your model.
year = 2005;
for mmm = 1:length(mloc_list)
    mloc = mloc_list{mmm};
    [mlon(mmm),mlat(mmm)] = Z_mooringLocations(mloc,year);
end
% The "tag" is just appended to the outfile name.
tag = 'RNtest';
% &&&&&&&&&&& END USER INPUT &&&&&&&&&&&&&&&

outfile = [Tdir.moor_out,basename,'_',tag,'.mat'];

% screen output
disp(' ')
disp('***************************************************')
disp(['indir = ',indir])
disp(['basename = ',basename])
for mmm = 1:length(mloc_list)
    disp(['mloc = ',mloc_list{mmm},', year = ',num2str(year)])
    disp([' ** lon = ',num2str(mlon(mmm)), ...
        ' lat = ',num2str(mlat(mmm))])
end
disp(['outfile = ',outfile])
disp('NOTE: run structure only saved in M(1)')
disp('---------------------------------------------------')

% get the run definition
run = roms_createRunDef('my_run',indir);

% extract the mooring
[M] = Z_get_moor(run,mlon,mlat);

% add to the results structure
M(1).run = run;
M(1).basename = basename;
for mmm = 1:length(mloc_list)
    mloc = mloc_list{mmm};
    M(mmm).mloc = mloc;
end

% save the results
save(outfile,'M');
disp('---------------------------------------------------')

