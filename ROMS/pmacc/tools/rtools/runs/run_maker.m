% run_maker.m 12/20/2012 Parker MacCready & Sarah Giddings
%
% driver for the code that creates ROMS forcing files
%
% Uses object-oriented programming.

clear all

% set global paths and directories
addpath('../../alpha/'); Tdir = toolstart;

% set local paths required by some methods
addpath run_definitions run_functions

% add user-specific paths
addpath([Tdir.tools 'rtools_user/runs/run_definitions']);


disp('**********************************')

gname = 'scidac'; % defines the grid to use
tag1 = '2'; % with the grid, defines the run definition to use
tag2 = '2005'; % an additional tag to denote the time interval
% and this is the actual time interval
%tdlims = [datenum(2005,1,1,0,0,0),datenum(2006,1,1,0,0,0)];
tdlims = [datenum(2005,12,24,0,0,0),datenum(2005,12,31,0,0,0)]; % TESTING

% initial creation of the run object
eval(['rn = ',gname,'_',tag1,';']);

rn.addInfo(Tdir,gname,tag1,tag2,tdlims);

rn.makeDir;

rn.addGrid;

rn.setScoord;

rn.makeS;

rn.makeClim;

rn.makeAtm;

rn.makeTide;

rn.makeRivers;

rn.addDye;
