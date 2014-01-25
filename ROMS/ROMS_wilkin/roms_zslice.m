function [data,x,y,t,grd] = roms_zslice(file,var,time,depth,grd)
% $Id: roms_zslice.m 422 2014-01-13 17:13:14Z wilkin $
% Get a constant-z slice out of a ROMS history, averages or restart file
% [data,x,y] = roms_zslice(file,var,time,depth,grd)
%
% Inputs
%    file = his or avg nc file
%    var = variable name
%    time = time index in nc file
%    depth = depth in metres of the required slice
%    grd (optional) is the structure of grid coordinates from roms_get_grid 
%
% Outputs
%    
%    data = the 2d slice at requested depth 
%    x,y = horizontal coordinates
%    t = time in days for the data
%
% John Wilkin

depth = -abs(depth);

% open the history or averages file
%nc = netcdf(file);

if ~nc_isvar(file,var)
  error([ 'Variable ' var ' is not present in file ' file])
end

% get the time

% this change (2013-05-23) to accommodate Forecast Model Run 
% Collection (FMRC) which changes the time coordinate to be named
% "time" but leaves the attribute of the variable pointed to ocean_time
info = nc_vinfo(file,var);
time_variable = info.Dimensions(end).Name;
% time_variable = nc_attget(file,var,'time');

if nc_varsize(file,time_variable)<time
  disp(['Requested time index ' int2str(time) ' not available'])
  disp(['There are ' int2str(nc_varsize(file,time_variable)) ...
    ' time records in ' file])
  error(' ')
end
t = roms_get_date(file,time); % gets output in matlab datenum convention

% check the grid information
if nargin<5 || (nargin==5 && isempty(grd))
  % no grd input given so try to get grd_file name from the history file
  grd_file = file;
  grd = roms_get_grid(grd_file,file);
else
  if ischar(grd)
    grd = roms_get_grid(grd,file);
  else
    % input was a grd structure but check that it includes the z values    
    if ~isfield(grd,'z_r')
      error('grd does not contain z values');
    end
  end
end

% get the data to be zsliced
%data = nc_varget(file,var,[time-1 0 0 0],[1 -1 -1 -1]);
data = ncread(file,var,[1 1 1 time],[Inf Inf Inf 1]);
data = permute(data,[3 2 1]);

% THIS STEP TO ACCOMMODATE NC_VARGET RETURNING A TIME LEVEL WITH
% LEADING SINGLETON DIMENSION - BEHAVIOR THAT DIFFERS BETWEEN JAVA AND
% MATLAB OPENDAP INTERFACES - 11 Dec, 2012
data = squeeze(data);

% slice at requested depth
[data,x,y] = roms_zslice_var(data,1,depth,grd);

switch roms_cgridpos(size(data),grd)
  case 'u'
    mask = grd.mask_u;
  case 'v'
    mask = grd.mask_v;
  case 'psi'
    mask = grd.mask_psi;
  case 'rho'
    mask = grd.mask_rho;
end

% Apply mask to catch shallow water values where the z interpolation does
% not create NaNs in the data

mask(mask == 0) = NaN;
data = data.*mask;
