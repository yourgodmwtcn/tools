% mooring_extractor.m  3/28/2012  Parker MacCready
%
% this is for extracting the results of ROMS simulations
% it is focused on time series at a point, over the full water column, as
% would be collected by a mooring
%
% uses linear interpolation to get the values at the right location
%
% NOTE: this calls a new version of the extraction code:
% [M] = Z_get_moor(run,mlon,mlat);
% which on my machine is much faster than the old version.

clear;
% &&&&&&&&&&& USER EDIT THIS &&&&&&&&&&&&&&
%indir = '/Users/PM/Documents/Salish/runs/';
indir = 'E:\Work\CattlePass\runs\';
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

[Tdir] = pan_start; % locate things and set paths
% create output directory if needed
moordir = [Tdir.pan_results,'moorings/'];
if ~exist(moordir,'dir'); mkdir(moordir); end;

% choose run to work on and set the basename
disp('********* mooring_extractor.m *******************')
disp(' ')
if 0 % interactive version
    disp('** Select OUT where history files are **');
    pth = uigetdir(indir);
else % hardwired version
    pth = [indir,'salish_2006_4/OUT'];
end
% assumes that "pth" is something like:
% /Users/PM/Documents/Salish/runs/ptx_med_2005_1/OUT
% then the lines below find the string right before "OUT", which is the
% "basename," the main identifier of the run
ind = strfind(pth,'/');
basename = pth(ind(end-1)+1:ind(end)-1);
disp(' '); disp(['basename = ',basename]); disp(' ')

% get the run definition
run = roms_createRunDef('my_run',pth);
% get the year
year = datestr(run.his.nctime(1),'yyyy');
year = str2num(year);

if 0 % choose a RISE mooring location
    mloc = 'rn'; % rs,rn,rc all valid
    % only good for 2004:2006
    [mlon,mlat] = Z_RISE_moor_loc(year,mloc);
else % USER supplied values
    % SciDAC
    mlon = -125.9;
    mlat = 48.5;
    mloc = 'scidac1';
    % Cha'ba/NEMO mooring
%     mlon = -125;
%     mlat = 48;
%     mloc = 'nemo_new';
    % Admiralty Inlet for Ocn Sci 2012 talk
    % mlon = -122.7;
    % mlon = 48.13;
    % mloc = 'ai0';
end

[M] = Z_get_moor(run,mlon,mlat);

% pack the results in a structure
M.run = run;
M.basename = basename;

% rename for backward consistency
mod_moor = M;

save([moordir,basename,'_',mloc,'.mat'],'mod_moor');
