% moor_start.m  8/20/2013  Parker MacCready
%
% sets up paths and directories (in structure Tdir) for mooring code

addpath('../alpha');
Tdir = toolstart;
Tdir.moor_out = [Tdir.output,'moor_out/'];
if ~exist(Tdir.moor_out,'dir'); mkdir(Tdir.moor_out); end;

% add path of user plot code
addpath('../moor_user/plot_code/');

addpath('Z_functions');
addpath('Z_plot_code');
