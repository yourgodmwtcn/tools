% Sets the appropriate lat/lon based on the year and mooring location
% all RISE and ECOHAB moorings are included
% as are Rick Thompson's J2C, A1, and DFO-MPO racerocks
% note all locations in 2004, 2005, and 2006 are included
% if a mooring was not deployed in one of those years, the location
% defaults to the 2006 location
% alternatively, a user can specify a lat/lon location and a mooring name
% other than those included here
% SNG 26 April 2012, updated 4 June 2012
%
% PM 6/8/2012 fixed RS 2004 lat_moor location
% also renamed by PM 6/8/2012 to fit into pandora
% NOTE: this should be recoded to be a function instead of a program

disp(['************ Year = ',num2str(yeari),' ***********']);
disp(['*** Mooring location: ' mloc])

switch mloc
    case 'RS'   % RISE South mooring
        switch yeari
            case 2004
                lat_moor = 46.0530; lon_moor = -124.1005;
                % fixed typo in lat_moor (was 45.0530) PM 6/8/2012
            otherwise %2005, 2006 locations
                lat_moor = 45.5; lon_moor = -124.1029;
        end
    case 'RC'   % RISE Central (plume) mooring
        switch yeari
            case 2004
                lat_moor = 46.1671; lon_moor = -124.1954;
            case 2005
                lat_moor = 46.1666; lon_moor = -124.1954;
            otherwise %2006 location
                lat_moor = 46.1668; lon_moor = -124.1955;
        end
    case 'RN'   % RISE North mooring
        switch yeari
            case 2004
                lat_moor = 46.4374; lon_moor = -124.3013;
            case 2005
                lat_moor = 46.9997; lon_moor = -124.4919;
            otherwise %2006 location
                lat_moor = 47.0166; lon_moor = -124.4921;
        end
    case 'EH1'  % ECOHAB JdF mouth mooring
        switch yeari
            case 2004
                lat_moor = 48 + (29.31/60); lon_moor = -(124 + (41.98/60));
            case 2005
                lat_moor = 48 + (29.330/60); lon_moor = -(124 + (41.956/60));
            otherwise %2006 location
                lat_moor = 48 + (29.391/60); lon_moor = -(124 + (41.919/60));
        end
    case 'EH2'  % ECOHAB south mooring
        switch yeari
            case 2004 % ECOHAB Coast Strait mooring (diff location from other years)
                lat_moor = 47 + (36.02/60); lon_moor = -(124 + (46.06/60));
            case 2005
                lat_moor = 47 + (35.789/60); lon_moor = -(124 + (35.906/60));
            otherwise %2006 location
                lat_moor = 47 + (35.861/60); lon_moor = -(124 + (35.415/60));
        end
    case 'EH4'  % another ECOHAB mooring inshore of EH2, not deployed in 2004, only 2005, 2006
        switch yeari
            case 2005
                lat_moor = 47 + (36.041/60); lon_moor = -(124 + (32.009/60));
            otherwise %2006 location
                lat_moor = 47 + (36.068/60); lon_moor = -(124 + (32.107/60));
        end
    case 'EH3'  % ECOHAB Eddy mooring
        switch yeari
            case 2004
                lat_moor = 48 + (17.82/60); lon_moor = -(125 + (27.54/60));
            case 2005
                lat_moor = 48 + (17.771/60); lon_moor = -(125 + (27.280/60));
            otherwise %2006 location
                lat_moor = 48 + (17.725/60); lon_moor = -(125 + (27.296/60));
        end
    case 'J2C'  % Rick Thomson (2007) JF2C Juan de Fuca center mooring, mooring moves around slightly, this is an average location
        lat_moor = 48 + (21.250/60); lon_moor = -(124 + (12.296/60));
    case 'A1'   % Rick Thomson IOS A1 offshore mooring (2005 location)
        lat_moor = 48.5293; lon_moor = -126.2025;
    case 'racerocks'    %data from http://www.pac.dfo-mpo.gc.ca/science/oceans/data-donnees/lighthouses-phares/index-eng.htm
        lat_moor = 48 + (17/60) + (54/(60*60)); lon_moor = -(123 + (31/60) + (54/(60*60)));
    otherwise
        if ~exist('lati','var')
            error([mloc ' is not defined in mooringLocations.m, you must specify the lat/lon location'])
        else
            lat_moor = lati; lon_moor = loni;
        end
end
disp(['*** lat = ' num2str(lat_moor) ', lon = ' num2str(lon_moor)])
