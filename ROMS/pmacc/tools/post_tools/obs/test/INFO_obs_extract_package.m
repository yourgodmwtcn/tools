%% examples and information for  OBS_EXTRACT package 
%{
CONTENTS:  ( see header info in .m files for more information)

obs_extract - extracts one dimensional (CTD or Bottle) data from multiple 
   files or folders and merges data  into one structure, will ignore 
    2-dimensional data.  Also omits data with nans in the coordinate
    varaibales.  Also uses obs_identifyCast to append 
   .cast to data structure with a unique number for each x-y-t triplet
obs_extractFromFile -  extracts one dimensional (CTD or Bottle) data from 
    a single file, used by obs_extract
obs_mooringExtract - extracts two dimensional data (Moorinig) data from a 
    single file,
obs_identifyCasts -finds unique x-y-t triplet in 1 dimensional data and 
    asigns a unique number,1-D data only
obs_seperateCasts - seperates data from a given cast into multiple seperate structures
    1-D data only
obs_merge - merges multiple data structures into one structure, 1-D data only
obs_omit - find and deletes data in structure  from all fields  where criterion is met
    works with both 1-D and 2-D
utility/trackDist - calculates distance in KM along a trackline as well as 
    closest point on a trackline to point off the line
utility/dirwalker - searches through specified folder for all netcdf files nested within folder


 HCthalweg.mat - Hood Canal Thalweg
    

.m files from other institutions included/needed by obs_extract package
    make_range_ring
    sw_pres
%}


% C. Bassin , cbassin@apl.washington.edu
% N. Banas and D. Sutherland  UW
% Mar 2010


%%  
%{
data = obs_extract(filename, vars, timeRange
                                        ..., 'section',x,y,range);
                                        ..., 'polygon',x,y);
                                        ..., 'all');
data = obs_extract(directory, ...
data = obs_extract(cell array of files and directories, ...


vars can be a string containing one variable name, or a cell array of strings,
  or 'all'.

timeRange is a 2 element vector of matlab datenums or 'all'.

'section'  use if interested in data within ? km from a thalweg or single
  location, x is vector of longitudes, y is vector of latitudes, range is
  distance in km from point or thalweg of interest
'polygon'  use if interested in data wihtin polygon,  x is vector of
  longitudes, y is vector of latitudes.
'all'  data from all locations returned

%}



%% flourescence data from the last decade within 1 km of Admirality Inlet station DOE ADM002 -122.6167,48.03
% when interested in data within proximity to one station, use 'section' option
files='..\MatlabWork\PSVS\';  % folder you want to search 

data = obs_extract(files, 'fluorescence', [datenum(2000,1,1) datenum(2010,1,1)], 'section',-122.6167,48.03,1);

%% CTD data for Hood Canal during Summer 2005 within 2 km of thalweg
files='..\MatlabWork\PSVS\' ; % folder you want to search 
load('thalweg.mat'); %load thalweg longitude and latitude data
variables_wanted={'temperature','salinity','fluorescence','density','oxygen'};
data = obs_extract(files,variables_wanted , [datenum(2005,6,1) datenum(2005,9,1)], 'section',LonVi,LatVi,2);




%% MOORING DATA EXAMPLE
% inputs are the same as obs_extract, except filename must be a string
% containing one file only.  Does not work on folders or mulitple files

filename='..\MatlabWork\ORCA_netcdf\ORCA_Twanoh.nc';
variables_wanted={'temperature','salinity','density','fluorescence','nitrate','oxygen'};
EX2D = obs_mooringExtract(filename, variables_wanted, [datenum(2006,1,1) datenum(2007,1,1)],'all');


