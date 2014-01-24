% lp_driver.m 3/24/2011 Parker MacCready
%
% This is code for making tidally-averaged versions of the history files
% from a ROMS run, using Neil's psvs code

% Notes from Neil:
% This will create a new series of files with suffix _lp.nc. You can use
% run.hislp the same way you use run.his in my other roms_extract
% examples.
% 
% This syntax might make more sense if you think of it in object terms.
% "run" is like a Run object, that contains FileSeries objects like
% run.his, run.hislp along with the grid, a nickname, and other
% metadata. If we were using the diagnostic files, there would be a
% run.dia. Later we might standardize a way of linking in forcing time
% series (run.forcing.wind, etc) while maintaining back-compatibility
% with code we write using roms_extract now.

% NOTE: I have edited the psvs/roms code
%   roms_hanning.m
%   roms_averageFileSet.m
% so that they accept an extra input argument "outdir" and all the output
% is written to this (presumably local) directory

% SNG moved user choices to beginning. I also adjusted runpath and rundir 
% so that you do not need to call this code from the run directory itself
% finally I added paths to run totally remotely on skua
% 6/20/2011 SNG update to include roms_hanning.m option of istart
% 4/16/2012 SNG updated roms_godinfilt.m to include istart and outdir and
% this file to be able to directy call roms_godinfilt
% 4/17/2012 SNG updated to use roms_godinfilt_native (using MATLAB's native
% netcdf language)

%% USER INPUT SECTION
fasttest = 0; %set to 0 to run a full lp, 1 to do a faster test

%set global paths and directories
%specifically make sure the path to your tools dir is set up within toolstart
addpath('../../alpha/'); Tdir = toolstart;

%set the runname
runname = 'ptx_highT_2_2005_addxxobcTout3602';
%set the run path
runpath = ['/pmraid3/sarahgid/runs/', runname];
%set the desired outpath (could be the same as the run path or different)
outpath = ['/pmraid3/sarahgid/runs/', runname];
%chose the output time base, the default [] is 1 file per day at noon
%alternatively you could input a vector as shown below
%outputTimebase = [];
outputTimebase = datenum('Jan 1 2005 12:00') : 1 : datenum('Dec 31 2005 12:00');
%set the start output file #, only change this if you wish to start mid-way through a run
istart = 1; 

%set the output directory to the outpath/OUT_lp/
%note the user could change suffix to change the output directory name if desired
suffix = 'lp';

%% set up the output directory and get run information
rundir = [runpath,'/OUT/'];
outdir = [outpath,'/OUT_',suffix,'/'];
if exist(outdir,'dir')==7
    disp(['Using existing output directory: ' outdir])
else
    mkdir(outdir)
end

%get run information
run = roms_createRunDef(rundir);

%% perform the low pass filter
if fasttest % toggle this for fast testing or regular tidal averaging
    % short window NOTE THIS IS FOR TESTING ONLY!
    run.hislp = roms_hanning(run.his,3/24,[],suffix,outdir);
else
    tic    
    %standard 40h hanning window
    % run.hislp = roms_hanning(run.his,40/24,outputTimebase,suffix,outdir,istart);
    %godin filter (24, 24, 25 hr)
    % run.hislp = roms_godinfilt(run.his,outputTimebase,suffix,outdir,istart);
    %godin filter (24, 24, 25 hr) with native netcdf language
    run.hislp = roms_godinfilt_native(run.his,outputTimebase,suffix,outdir,istart);
    toc
end

