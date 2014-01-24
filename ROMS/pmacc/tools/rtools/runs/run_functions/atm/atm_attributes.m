function [ncvarname,nclongname,ncunits,nctimename, ...
    scalefactor,scalefactor2] = atm_attributes(var);
% atm_attributes.m  12/21/2012  Parker MacCready
%
% used by make_atm.m

scalefactor = 1; % multiply by this
scalefactor2 = 0; % and then add this
switch var
    case 'psfc'
        ncvarname = 'Pair';
        nclongname = 'surface air pressure';
        ncunits = 'millibar';
        nctimename = 'pair_time';
        scalefactor = 1/100; % convert Pa to mbar
    case 'rain'
        ncvarname = 'rain';
        nclongname = 'rain fall rate';
        ncunits = 'kilograms meter-2 second-2';
        nctimename = 'rain_time';
        scalefactor = 10; % convert cm s-1 to kg m-2 s-1
    case 'swdown'
        ncvarname = 'swrad';
        nclongname = 'solar shortwave radiation flux';
        ncunits = 'watts meter-2';
        nctimename = 'srf_time';
        scalefactor = 1 - 0.1446; % account for reflection
    case 'lwdown'
        ncvarname = 'lwrad_down';
        nclongname = 'downwelling longwave radiation flux';
        ncunits = 'watts meter-2';
        nctimename = 'lrf_time';
    case 't2'
        ncvarname = 'Tair';
        nclongname = 'surface air temperature';
        ncunits = 'Celsius';
        nctimename = 'tair_time';
        scalefactor2 = -273.15; % convert K to C
    case 'qair'
        ncvarname = 'Qair';
        nclongname = 'surface air relative humidity';
        ncunits = 'percentage';
        nctimename = 'qair_time';
    case 'u10r'
        ncvarname = 'Uwind';
        nclongname = 'surface u-wind component';
        ncunits = 'meter second-1';
        nctimename = 'wind_time';
    case 'v10r'
        ncvarname = 'Vwind';
        nclongname = 'surface v-wind component';
        ncunits = 'meter second-1';
        nctimename = 'wind_time';
end