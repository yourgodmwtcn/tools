function [coast_handle] = Z_addcoast(whichfile,coastdir)
% 5/9/2012  Parker MacCready
%
% [coast_handle] = Z_addcoast(whichfile,coastdir)
%
% This adds a coastline to an existing plot, also does river tracks as of
% 5/12/2011 (need a better data file for this!)
%
% whichfile is 'detailed', 'regional', 'combined', or 'rivers'
%
% typically coastdir would be Tdir.coast (from alpha/toolstart.m)

hold on
switch whichfile
    case 'detailed'
        load([coastdir,'pnw_coast_detailed.mat']);
    case 'regional'
        load([coastdir,'pnw_coast_regional.mat']);
    case 'combined'
        load([coastdir,'pnw_coast_combined.mat']);
    case 'rivers'
        load([coastdir,'pnw_rivers.mat']);
end
coast_handle = plot(lon,lat,'-k');
hold off

