function [Tdir] = toolstart
% toolstart.m  6/6/2012  Parker MacCready
%
% A function to be invoked at the start of any primary code in the "tools"
% system.
%
% It sets required paths and returns a structure "Tdir" with locations of
% useful directories.
%
% This is intended to be hopefully the ONLY place a user of the tools code
% collections would have to set machine-dependent paths and file locations
%
% copyright Parker MacCready 2011, released under the BSD license

% &&&&&&&&& USER EDIT THIS IF NEEDED &&&&&&&&&&&&&&&&
% note that this is designed so that it should not have to
% be changed at all when moving to different systems, as long
% as the suggested directory structure is maintained
this_dir = pwd; t_ind = strfind(this_dir,'\tools\pandora');
Tdir.tools_parent = this_dir(1:t_ind);
Tdir.tools_parent = 'E:\Work\tools\ROMS\pmacc\';
% NOTE 6/1/2012 The code as written above should work from anywhere within
% the directories tools, tools_data, tools_output, and anywhere alse whose
% name starts with tools - as long as it is at the same level.
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% these may not have to be changed at all if you are using the suggested
% tools directory structure
%
% upper level places:
Tdir.tools = [Tdir.tools_parent,'tools\'];
Tdir.data = [Tdir.tools_parent,'tools_data\'];
Tdir.output = [Tdir.tools_parent,'tools_output\'];
% make sure that the output directory exists
if ~exist(Tdir.output,'dir'); mkdir(Tdir.output); end;

% more specific places:
Tdir.rtools = [Tdir.tools,'rtools/'];
% data places
Tdir.coast = [Tdir.data,'mossea_geo_data/coast/'];
Tdir.topo = [Tdir.data,'mossea_geo_data/topo/'];
Tdir.atm = [Tdir.data,'mossea_forcing_data/atm/'];
Tdir.river = [Tdir.data,'mossea_forcing_data/river/'];
Tdir.tide = [Tdir.data,'mossea_forcing_data/tide/'];
Tdir.ocn = [Tdir.data,'mossea_forcing_data/ocn/'];

% paths to shared code assumed to be available by many programs
addpath([Tdir.rtools,'Z_utils']);
%addpath([Tdir.tools,'shared/mexcdf/mexnc']);
%addpath([Tdir.tools,'shared/mexcdf/snctools']);
%addpath([Tdir.tools,'shared/seawater']);
addpath([Tdir.tools,'pandora/Z_functions']);
addpath([Tdir.tools,'post_tools/obs']);
addpath([Tdir.tools,'post_tools/roms']);
addpath([Tdir.tools,'post_tools/utility']);


