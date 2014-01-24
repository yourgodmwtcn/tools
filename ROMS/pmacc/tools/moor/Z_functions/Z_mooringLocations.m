function [lon_moor,lat_moor] = Z_mooringLocations(mloc,yeari)
%
% Sets the appropriate lat/lon based on the year and mooring location
% all RISE and ECOHAB moorings are included, as are Rick Thompson's
% J2C, A1, and DFO-MPO racerocks.  Also more locations as detailed in the
% change log below.
%
% ** If a mooring was not deployed in one of a given years, the location
% sometimes defaults to the 2006 location. **
%
% SNG 26 April 2012, updated 4 June 2012
% PM 6/8/2012 fixed RS 2004 lat_moor location
% PM 11/20/2012 Recoded as a function
% PM 1/13/_2013 Renamed as a Z function
% PM 2/8/2013 Added many OCNMS locations, and deepshelf47, BCshelf,
%    and ORshelf, from Sam Siedlecki
%
% INPUT:
%   mloc is a string, e.g. 'RS' for RISE South
%   yeari is the year, e.g. 2005 (a float)
%
% OUTPUT:
%   lon_moor (degrees) mooring longitude, -180 to 180 format
%   lat_moor (degrees) mooring latitude

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
                lat_moor = 48 + (29.31/60);
                lon_moor = -(124 + (41.98/60));
            case 2005
                lat_moor = 48 + (29.330/60);
                lon_moor = -(124 + (41.956/60));
            otherwise %2006 location
                lat_moor = 48 + (29.391/60);
                lon_moor = -(124 + (41.919/60));
        end
    case 'EH2'  % ECOHAB south mooring
        switch yeari
            case 2004 % ECOHAB Coast Strait mooring
                % (different location from other years)
                lat_moor = 47 + (36.02/60);
                lon_moor = -(124 + (46.06/60));
            case 2005
                lat_moor = 47 + (35.789/60);
                lon_moor = -(124 + (35.906/60));
            otherwise %2006 location
                lat_moor = 47 + (35.861/60);
                lon_moor = -(124 + (35.415/60));
        end
    case 'EH4'  % another ECOHAB mooring inshore of EH2,
        % not deployed in 2004, only 2005, 2006
        switch yeari
            case 2005
                lat_moor = 47 + (36.041/60);
                lon_moor = -(124 + (32.009/60));
            otherwise %2006 location
                lat_moor = 47 + (36.068/60);
                lon_moor = -(124 + (32.107/60));
        end
    case 'EH3'  % ECOHAB Eddy mooring
        switch yeari
            case 2004
                lat_moor = 48 + (17.82/60);
                lon_moor = -(125 + (27.54/60));
            case 2005
                lat_moor = 48 + (17.771/60);
                lon_moor = -(125 + (27.280/60));
            otherwise %2006 location
                lat_moor = 48 + (17.725/60);
                lon_moor = -(125 + (27.296/60));
        end
    case 'J2C'  % Rick Thomson (2007) JF2C Juan de Fuca center mooring,
        % mooring moves around slightly, this is an average location
        lat_moor = 48 + (21.250/60);
        lon_moor = -(124 + (12.296/60));
    case 'A1'   % Rick Thomson IOS A1 offshore mooring (2005 location)
        lat_moor = 48.5293; lon_moor = -126.2025;
    case 'racerocks'
        % data from http://www.pac.dfo-mpo.gc.ca/science/oceans/ ...
        % data-donnees/lighthouses-phares/index-eng.htm
        lat_moor = 48 + (17/60) + (54/(60*60));
        lon_moor = -(123 + (31/60) + (54/(60*60)));
        
        % ** START new listings from Samantha Siedlecki 2/8/2013 **
    case 'deepshelf47'   % ECOHAB SS mooring
        lat_moor = 47 + (36.041/60); lon_moor = -(125);
    case 'BCshelf'
        lat_moor= 49; lon_moor=-125.8;
    case 'ORshelf'
        lat_moor=45; lon_moor=-124.1;
    case 'CA015'
        switch yeari
            case 2009
                lat_moor=48+(9/60)+(58.68/(60*60));
                lon_moor=-(124+(45/60)+(24.6/(60*60)));
            case 2006
                lat_moor=48+(9/60)+(54/(60*60));
                lon_moor=-(124+(45/60)+(24/(60*60)));
            case 2007
                lat_moor=48+(9/60)+(58.68/(60*60));
                lon_moor=-(124+(45/60)+(24.6/(60*60)));
        end
    case 'CA042'
        switch yeari
            case 2009
                lat_moor=48+(9/60)+(58.08/(60*60));
                lon_moor=-(124+(48/60)+(14.58/(60*60)));
            case 2006
                lat_moor=48+(9/60)+(54/(60*60));
                lon_moor=-(124+(48/60)+(12/(60*60)));
            case 2007
                lat_moor=48+(9/60)+(54/(60*60));
                lon_moor=-(124+(48/60)+(12/(60*60)));
        end
    case 'CA065'
        switch yeari
            case 2009
                lat_moor=48+(9/60)+(57.12/(60*60));
                lon_moor=-(124+(53/60)+(41.64/(60*60)));
            case 2006
                lat_moor=48+(9/60)+(54/(60*60));
                lon_moor=-(124+(53/60)+(35/(60*60)));
            case 2007
                lat_moor=48+(9/60)+(54/(60*60));
                lon_moor=-(124+(53/60)+(34.8/(60*60)));
        end
    case 'CA100'
        lat_moor=48+(9/60)+(57/(60*60));
        lon_moor=-(124+(55/60)+(55/(60*60)));
    case 'CE015'
        switch yeari
            case 2009
                lat_moor=47+(21/60)+(15/(60*60));
                lon_moor=-(124+(20/60)+(56.58/(60*60)));
            case 2006
                lat_moor=47+(21/60)+(19/(60*60));
                lon_moor=-(124+(21/60)+(8/(60*60)));
            case 2007
                lat_moor=47+(21/60)+(19.2/(60*60));
                lon_moor=-(124+(21/60)+(7.8/(60*60)));
        end
    case 'CE042'
        switch yeari
            case 2009
                lat_moor=47+(21/60)+(11.28/(60*60));
                lon_moor=-(124+(29/60)+(19.44/(60*60)));
            case 2006
                lat_moor=47+(21/60)+(19/(60*60));
                lon_moor=-(124+(28/60)+(6/(60*60)));
            case 2007
                lat_moor=47+(21/60)+(19.2/(60*60));
                lon_moor=-(124+(28/60)+(6/(60*60)));
        end
    case 'CE065'
        switch yeari
            case 2009
                lat_moor=47+(21/60)+(10.2/(60*60));
                lon_moor=-(124+(34/60)+(0.96/(60*60)));
            case 2006
                lat_moor=47+(21/60)+(19/(60*60));
                lon_moor=-(124+(34/60)+(7/(60*60)));
            case 2007
                lat_moor=47+(21/60)+(19.2/(60*60));
                lon_moor=-(124+(34/60)+(7.2/(60*60)));
        end
    case 'CE100'
        lat_moor=47+(21/60)+(19/(60*60));
        lon_moor=-(124+(40/60)+(34/(60*60)));
    case 'KL015'
        switch yeari
            case 2009
                lat_moor=47+(36/60)+(3/(60*60));
                lon_moor=-(124+(25/60)+(42.24/(60*60)));
            case 2007
                lat_moor=47+(35/60)+(54/(60*60));
                lon_moor=-(124+(25/60)+(36/(60*60)));
            case 2006
                lat_moor=47+(35/60)+(54/(60*60));
                lon_moor=-(124+(25/60)+(36/(60*60)));
        end
    case 'KL027'
        switch yeari
            case 2009
                lat_moor=47+(35/60)+(40.44/(60*60));
                lon_moor=-(124+(29/60)+(49.44/(60*60)));
            case 2007
                lat_moor=47+(35/60)+(36/(60*60));
                lon_moor=-(124+(29/60)+(36/(60*60)));
            case 2006
                lat_moor=47+(35/60)+(36/(60*60));
                lon_moor=-(124+(29/60)+(36/(60*60)));
        end
    case 'MB015'
        switch yeari
            case 2009
                lat_moor=48+(19/60)+(31.44/(60*60));
                lon_moor=-(124+(40/60)+(54.54/(60*60)));
            case 2007
                lat_moor=48+(19/60)+(30/(60*60));
                lon_moor=-(124+(40/60)+(54/(60*60)));
            case 2006
                lat_moor=48+(19/60)+(30/(60*60));
                lon_moor=-(124+(40/60)+(54/(60*60)));
        end
    case 'MB042'
        switch yeari
            case 2009
                lat_moor=48+(19/60)+(26.28/(60*60));
                lon_moor=-(124+(44/60)+(7.38/(60*60)));
            case 2007
                lat_moor=48+(19/60)+(18/(60*60));
                lon_moor=-(124+(44/60)+(36/(60*60)));
            case 2006
                lat_moor=48+(19/60)+(18/(60*60));
                lon_moor=-(124+(44/60)+(36/(60*60)));
        end
    case 'TH015'
        switch yeari
            case 2009
                lat_moor=47+(52/60)+(31.8/(60*60));
                lon_moor=-(124+(37/60)+(7.14/(60*60)));
            case 2007
                lat_moor=47+(52/60)+(30/(60*60));
                lon_moor=-(124+(37/60)+(6/(60*60)));
            case 2006
                lat_moor=47+(52/60)+(30/(60*60));
                lon_moor=-(124+(37/60)+(6/(60*60)));
        end
    case 'TH042'
        switch yeari
            case 2009
                lat_moor=47+(52/60)+(34.14/(60*60));
                lon_moor=-(124+(44/60)+(0.3/(60*60)));
            case 2007
                lat_moor=47+(52/60)+(34.2/(60*60));
                lon_moor=-(124+(44/60)+(0.6/(60*60)));
            case 2006
                lat_moor=47+(52/60)+(34/(60*60));
                lon_moor=-(124+(44/60)+(0/(60*60)));
        end
        % ** END new listings from Samantha Siedlecki 2/8/2013 **
        
    otherwise
        error([mloc ' is not defined in mooringLocations'])
end
