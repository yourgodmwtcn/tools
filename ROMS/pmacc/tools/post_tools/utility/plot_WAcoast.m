function plot_WAcoast(type, varargin)
%-----------------------------
%
% plot_WAcoast(type)
%
% this adds a coastline to the current plot of type:
% 'detailed': loads coastline_detailed.mat which is good for Puget Sound
% 'regional': loads coastline_regional.mat which is good for WA coast and S of Georgia
% 'regional_wbathy': includes regional coastline and contours at 50, 100,
% 200, 500, and 1000m depth
%
% varargin are line specs to pass to plot 
%
% ex) 'linewidth',2,'color','k' (default is black)
% das 5/18/2009
% sng updated 4/12/2011 to include bathymetry contours
%-----------------------------

if nargin < 1 %default to regional
    type = 'regional';
end

% define where these coastline files exist
addpath(genpath('/Users/sarahgid/Documents/romstools/mossea_geo_data/'));
dirname = '/Users/sarahgid/Documents/romstools/mossea_geo_data/';
%dirname = '/skua1/sarahgid/tools/mossea_geo_data/';
switch type
    case 'detailed'
        filename = [dirname, 'coastlines/coastline_detailed.mat'];
        load(filename);
    case 'regional'
        filename = [dirname, 'coastlines/coastline_regional.mat'];
        load(filename);
    case 'all'
        filename = [dirname, 'coastlines/coastline_combined.mat'];
        load(filename);
        lon_coast = lon;lat_coast = lat;
    case 'all_wbathy'
        filename = [dirname, 'coastlines/coastline_combined.mat'];
        load(filename);
        lon_coast = lon;lat_coast = lat;
        filename = [dirname, 'hybrid/pnw_combined_full.mat'];
        load(filename);
        contour(lon_topo,lat_topo,z_topo,[-4000 -3000 -2000 -1000 -500 -180 -100 -50 -30],'color',[0.5 0.5 0.5])
    case 'regional_wbathy'
        filename = [dirname, 'coastlines/coastline_regional.mat'];
        load(filename);
        filename = [dirname, 'hybrid/pnw_combined_full.mat'];
        load(filename);
        contour(lon_topo,lat_topo,z_topo,[-4000 -3000 -2000 -1000 -500 -180 -100 -50 -30],'color',[0.5 0.5 0.5])
    otherwise
        error('specify either detailed, regional, regional_wbathy, or nothing')
end

hold on;
plot(lon_coast, lat_coast, 'color', 'k', varargin{:})


