function [lon_moor,lat_moor] = Z_RISE_moor_loc(year,mloc)
% 12/14/2011  Parker MacCready
%
% gives RISE mooring locations
switch year
    case 2004
        switch mloc
            case 'rs'    % South mooring
                lat_moor = 46 + 03.179/60; lon_moor = -(124 + 06.027/60);
            case 'rc'    % Central (plume) mooring
                lat_moor = 46 + 10.024/60; lon_moor = -(124 + 11.724/60);
            case 'rn'    % North mooring
                lat_moor = 46 + 26.245/60; lon_moor = -(124 + 18.080/60);
        end
    case 2005
        switch mloc
            case 'rs'    % South mooring
                lat_moor = 45.5; lon_moor = -124.1029;
            case 'rc'    % Central (plume) mooring
                lat_moor = 46.1666; lon_moor = -124.1954;
            case 'rn'    % North mooring
                lat_moor = 46.9997; lon_moor = -124.4919;
        end
    case 2006
        switch mloc
            case 'rs'    % South mooring  (90 m)
                lat_moor = 45 + 30.002/60; lon_moor = -(124 + 06.176/60);
            case 'rc'    % Central (plume) mooring  (71 m)
                lat_moor = 46 + 10.007/60; lon_moor = -(124 + 11.731/60);
            case 'rn'    % North mooring (70 m)
                lat_moor = 47 + 00.995/60; lon_moor = -(124 + 29.525/60);
        end
end