% grid_maker.m 11/28/2012 Parker MacCready, Sarah Giddings, & Neil Banas
%
% driver for the code that creates a ROMS grid file
%
% Uses object-oriented programming.

clear all

disp('*******************')
% *** USER: Set the name of the grid here ***
gname = 's2';
% *******************************************

% set general paths
addpath('../../alpha/');
Tdir = toolstart;

% add paths specific to making grids
addpath([Tdir.rtools 'grid/grid_definitions']);
addpath([Tdir.rtools 'grid/grid_functions']);
addpath([Tdir.rtools 'grid/grid_functions/mask']);
addpath([Tdir.rtools 'grid/grid_functions/LP_Bathymetry/Mfiles']);

% add user-specific paths
addpath([Tdir.tools 'rtools_user/grid/grid_definitions']);

% initial creation of the grid object
eval(['grd = ',gname,';']);

% execute available methods as desired

% adding info the the object
grd.addGname(gname);
grd.addTdir(Tdir);
grd.addRiverFile;
grd.addTopoFiles;

% make the initial bathymetry and etc.
grd.makeLonLat;
grd.makeHraw;
grd.makeGridExtras;

% mask and carve rivers
grd.addMaskFile;
if ~isempty(grd.maskFile) % use an existing mask
    grd.carveRivers;
    grd.loadMask;
else % make and edit a new mask
    grd.makeMask;
    grd.carveRivers;
    grd.editMask;
end

% smooth and write to NetCDF
grd.smooth;
grd.writeNetCDF;

% plot results
plot_mask(grd,'h',Tdir)



