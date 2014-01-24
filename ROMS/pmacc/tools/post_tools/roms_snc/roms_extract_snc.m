function [data,coords] = roms_extract_snc(varargin)

% [data,coords] = roms_extract_snc(series, 3Dvarname, t, 'full');
%                                               ..., 'point', y, x);
%
% [data,coords] = roms_extract_snc(series, 4Dvarname, t, 'full');
%                                               ..., 'surface');
%                                               ..., 'bottom');
%                                               ..., 'zslice', z);
%                                               ..., 'depthslice', depth);
%                                               ..., 'depthaverage', [mindepth maxdepth]);
%                                               ..., 'depthintegral', [mindepth maxdepth]);
%                                               ..., 'profile', y, x);
%                                               ..., 'point', z, y, x);
%
%           ... = roms_extract_snc(filename, varname, ...
%
% general routine for extracting data from a ROMS netcdf file series or a
% single file (e.g., ocean_his_*.nc) in data units.
%
% *** this version uses the snctools rather than native matlab netcdf support. Since this
% makes the dimensions come out in a different order, it's a major change! The _snc*
% versions are just kept around for back compatibility (analysis code from mid 2013 or
% before). ***
%
% assumes that the first dimension is ocean_time, and that there's only one save
% per file, i.e.
%     size [1 J I] for 3D variables, where (J,I) match the rho, u, or v grid
%     size [1 K J I] for 4D variables, where K matches the rho or w grid, or K=1
%
% If a file series is specified (see roms_createSeriesDef.m), t can be
% either a scalar or a vector. If a filename is given, t isn't given at all.
% z, y, x can be scalars or vectors.
%
% zslice = measured from MSL; depthslice,average,integral = measured from surface.
%
% coords is a structure containing plaid matrices the same size as data:
%                  t scalar:       t vector:
%     3D var:      ym,xm           tm,ym,xm
%     4D var:      zm*,ym,xm       tm,zm,ym,xm
%
%                  *not for depthaverage,
%				   depthintegral
%
% neil banas mar 2009

if ischar(varargin{1})
	[data,coords] = roms_extractFromFile_snc(varargin{:});
elseif isstruct(varargin{1})
	[data,coords] = roms_extractFromSeries_snc(varargin{:});
end
