%
% Bathymetry Extracting Functions
% ===============================
%
% These functions are used for extracting bathymetry from a dataset. 
%
%
%   etopo4.nc     - ETOPO5 Earth surface topography dataset
%                   (5 minute resolution).
%
%   c_bath        - Create a bathymetry NetCDF file.
%   extract_bath  - Driver to extract bathymetry.  The extracted data is
%                   either written into a Matlab file that can be use in
%                   "seagrid" or a NetCDF file.
%   x_etopo5      - Extract requested bathymetry from ETOPO-5 dataset.
%

% svn $Id: Contents.m 586 2012-01-03 20:19:25Z arango $
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
