%
%EASYGRID   Generates grid for ROMS.
%
%   This is a "quick-and-dirty" way to generate a grid for ROMS, as
%   well as a simple initial conditions file where specified values
%   are constant over the grid domain.
%
%   USAGE: Before running this script, you need to specify the 
%          parameters in the "USER SETTINGS" section of this code. Beside
%          each variable there are detailed instructions. NOTE the 
%          "SWITCHES" section, which turns ON or OFF different
%          functionalities of EASYGRID, for example: edit_mask turns on/off 
%          the mask-editing routine and smooth_grid turns on/off bathymetry 
%          smoothing. You probably will need to run EASYGRID several times in 
%          a iterative manner to create, smooth and mask-edit a grid. See
%          tutorial link below.       
%
%   TUTORIAL: https://www.myroms.org/wiki/index.php/easygrid
%   
%   DEPENDENCIES: MEXNC (specific for your architecture), SNCTOOLS
%                 and the ROMS Matlab toolkit.
%                 For a tutorial on how to download/install dependencies,
%                 see: https://www.myroms.org/wiki/index.php/MEXNC
%
%   LIMITATIONS: This script should NOT be used to generate GLOBAL grids.
%
%   OUTPUTS: 1) NetCDF file (.nc) with grid for ROMS (if save_grid is ON).
%            2) NetCDF file (.nc) with initial conditions (if save_init is ON).
%            3) Screen output of many parameters that need to be 
%               copy-pasted into the ocean.in file (if screen_output is ON).
%            4) Plot of the grid (if plot_grid is ON).
%            5) Mask-Editing interactive plot (if edit_mask is ON), 
%               to change sea-cells into land-cells and vice versa.
%
%   TESTING: After download, run EASYGRID with all the default parameters.
%            If working properly, EASYGRID should save two NetCDF files 
%            (grid and initialization files) for the ROMS "Fjord Tidal Case".
%            Also screen-output of grid parameters and a plot of the grid
%            should be generated.
%
%==================================================================
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%==================================================================
%                                     Version: 0 Date: 2008-Apr-24
%
%                      Written by: Diego A. Ibarra (dibarra@dal.ca)
%             with borrowed code from Katja Fennel and John Wilkin,
%                 and with help of many functions by Hernan Arango.
% 
%   See also   SMTH_BATH   SPHERIC_DIST   PCOLORJW.


%****************************************************************************************************************************************************
% USER SETTINGS *************************************************************************************************************************************
%****************************************************************************************************************************************************

% -------------------------------------------------------------------------
% SWITCHES ( ON = 1; OFF = 0 ) --------------------------------------------
% -------------------------------------------------------------------------
    create_grid = 1;   % Create GRID. Turn OFF to work with a previously created grid (i.e. grid variables existing on Workspace)
      plot_grid = 1;   % Plot grid
    smooth_grid = 1;   % Smooth bathymetry using H. Arango's smth_bath.m , which applies a Shapiro filter to the bathymetry
      edit_mask = 0;   % Edit rho-mask using interactive plot. Use this to manually change sea-pixels into land-pixels and vice-versa 
  screen_output = 1;   % Displays (on screen) many parameters that need to be copy-pasted in the ocean.in file    
      save_grid = 1;   % Save GRID in a NetCDF file
      save_init = 1;   % Create (and save) INITIAL CONDITIONS (from grid)
                       % ON = 1; OFF = 0
    
% -------------------------------------------------------------------------
% GRID SETTINGS -----------------------------------------------------------
% -------------------------------------------------------------------------

% Geographical and Grid parameters --------
     lat =  44.8125;        % Latitude  (degrees) of the bottom-left corner of the grid.
     lon = -62.8855;        % Longitude (degrees) of the bottom-left corner of the grid. 

       X = 9100;            % Width of domain (meters)
       Y = 2550;            % Length of domain (meters)
rotangle = -43;             % Angle (degrees) to rotate the grid conterclock-wise
   resol = 100;              % Cell width and height (i.e. Resolution)in meters. Grid cells are forced to be (almost) square.
       N = 10;              % Number of vertical levels
     
       
% Bathymetry -------------- % Bathymetry for ROMS is measured positive downwards (zeros are not allowed) see: https://www.myroms.org/wiki/index.php/Grid_Generation#Bathymetry
                            % To allow variations in surface elevation (eg. tides) while keeping all depths positive,
                            % an arbitrary offset (see minh below) is added to the depth vector.
      
      hh = nan;             % Analytical Depth (meters) used to create a uniform-depth grid. If using a bathymetry file, leave hh = nan;
    minh = 4;               % Arbitrary depth offset in meters (see above). minh should be a little more than the maximum expected tidal variation.
   Bathy = 'Fjord_bathy.mat';% Bathymetry file. If using the analytical depth above (i.e. hh ~= nan), then Bathy will not be used.
                            % The bathymetry file should be a .mat file containing 3 vectors (xbathy, ybathy and zbathy). where xbathy = Longitude, 
                            % ybathy = Latitude and zbathy = depth (positive, in meters). xbathy and ybathy are in decimal degrees See: http://en.wikipedia.org/wiki/Decimal_degrees
       %Bathymetry smoothing routine...  See "Switches section" (above) to turn this ON or OFF
       if smooth_grid == 1;
           order = 2;       % Order of Shapiro filter (2,4,8)... default: 2
            rlim = 0.35;    % Maximum r-factor allowed (0.35)... default: 0.35
           npass = 50;      % Maximum number of passes.......... default: 50
       end
%---------------------------------
                       

% Coastline ----------------------
   Coast = load('Fjord_coast.mat'); % If there isn't a coastline file... comment-out this line: e.g. %Coast = load('Fjord_coast.mat');
                                    % The coastline is only used for plotting. The coastline .mat file should contain 2 vectors named "lat" and "lon"
                                     
       
% -------------------------------------------------------------------------
% OUTPUT: File naming and NetCDF descriptors ------------------------------
% -------------------------------------------------------------------------
Grid_filename = 'Fjord'; 	   % Filename for grid and initial conditions files (don't include extension). 
                               % "_grd.nc" is added to grid file and "_ini.nc" is added to initial conditions file
Descrip_grd   = 'Test grid';               %Description for grid .nc file
Descrip_ini   = 'Test initial conditions'; %Description for initial conditions .nc file
Author        = 'John Smith';
Computer      = 'My Computer';

% -------------------------------------------------------------------------
% INITIAL CONDITIONS ------------------------------------------------------
% -------------------------------------------------------------------------
if save_init == 1; % See "Switches section" (above) to turn this ON or OFF
    % Initial conditions will be constant throughout the grid domain
    %--------------------------------------------------------------------------
    NH4          = 0.1;     % Ammonium concentration (millimole_NH4 meter-3)
    NO3          = 10;      % Nitrate concentration (millimole_N03 meter-3)
    chlorophyll1 = 0.3;     % Chlorophyll concentration in small phytoplankyon (milligrams_chlorophyll meter-3)
    chlorophyll2 = 0.3;     % Chlorophyll concentration in large phytoplankyon (milligrams_chlorophyll meter-3)
    detritus1    = 0.03;    % Small detritus concentration (millimole_nitrogen meter-3)
    detritus2    = 0.03;    % Large detritus concentration (millimole_nitrogen meter-3)
    detritusC1   = 1;       % Small detritus carbon concentration (millimole_carbon meter-3)
    detritusC2   = 0.2;     % Large detritus carbon concentration (millimole_carbon meter-3)
    phyto1       = 0.02;    % Small phytoplankton concentration (millimole_nitrogen meter-3)
    phyto2       = 0.02;    % Large phytoplankton concentration (millimole_nitrogen meter-3)
    phytoC1      = 0.2;     % Small phytoplankton carbon concentration (millimole_carbon meter-3)
    phytoC2      = 0.1;     % Small phytoplankton carbon concentration (millimole_carbon meter-3)
    salt         = 30;      % Salinity (PSU)
    temp         = 9;       % Potential temperature (Celsius)
    u            = 0;       % u-momentum component (meter second-1)
    ubar         = 0;       % Vertically integrated u-momentum component (meter second-1)
    v            = 0;       % v-momentum component (meter second-1)
    vbar         = 0;       % Vertically integrated v-momentum component (meter second-1)
    zeta         = 0;       % Free-surface (meters)
    zooplankton  = 0.01;    % Zooplankton concentration "millimole_nitrogen meter-3"
    zooplanktonC = 0.5;     % Zooplankton carbon concentration "millimole_carbon meter-3"
    %--------------------------------------------------------------------------
end



%****************************************************************************************************************************************************
% END OF USER SETTINGS ******************************************************************************************************************************
%****************************************************************************************************************************************************


tic % Start the timer


disp([char(13),char(13),char(13),char(13),char(13),char(13)])
disp('***************************************************************')
disp('** EASYGRID ***************************************************')
disp('***************************************************************')


%Lets start with some conversions 
rotangle = rotangle/180*pi;                 % Convert Angle for grid rotation from degrees to radians
latdist  = spheriq_dist(lon,lat,lon,lat+1); % Length (in meters) of 1 degree of latitude




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START of GRID generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if create_grid == 1; % Only create grid if switch (in USER SETTINGS) is turned ON

    disp([char(13),'GENERATING grid...'])

    % Clear varibles from previous grids
    clear x y Lm Lp L Mm Mp M x_rho y_rho lat_rho lon_rho lat_u lon_u lat_v lon_v lat_psi lon_psi...
        el xl dx dy pm pn dndx dmde f h mask_rho mask_rho_nan mask_u mask_v mask_psi

    x = 0:resol:X; Lm = length(x)-2; Lp= Lm +2; L = Lp-1;
    y = 0:resol:Y; Mm = length(y)-2; Mp= Mm +2; M = Mp-1;

    % RHO GRID ------------------------------------------------------------
    %Create non-georeferenced grid in meters (origin = 0,0)
    x_rho = ones(length(y),1)*x(:)';
    y_rho = y(:)*ones(1,length(x));

    %Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_rho = x_rho * cos(rotangle) - y_rho * sin(rotangle);
    Ry_rho = x_rho * sin(rotangle) + y_rho * cos(rotangle);

    %Estimate Latitude and Longitude of each Grid point
    lat_rho(1,1) = lat + (Ry_rho(1,1) ./ latdist);
    lon_rho(1,1) = lon + (Rx_rho(1,1) ./ spheriq_dist(lon,lat_rho(1,1),lon+1,lat_rho(1,1)));
    for i = 2:length(Ry_rho(:,1));
        lat_rho(i,1) = lat + (Ry_rho(i,1) ./ latdist);
        lon_rho(i,1) = lon_rho(i-1,1) + ((Rx_rho(i,1)-Rx_rho(i-1,1)) ./ spheriq_dist(lon_rho(i-1,1),lat_rho(i,1),lon_rho(i-1,1)+1,lat_rho(i,1)));
    end

    for i = 1:length(Ry_rho(:,1));
        for j = 2:length(Ry_rho(1,:));
            lat_rho(i,j) = lat + (Ry_rho(i,j) ./ latdist);
            lon_rho(i,j) = lon_rho(i,j-1) + ((Rx_rho(i,j)-Rx_rho(i,j-1)) ./ spheriq_dist(lon_rho(i,j-1),lat_rho(i,j),lon_rho(i,j-1)+1,lat_rho(i,j)));
        end
    end
    clear Rx_rho Ry_rho

    % U GRID --------------------------------------------------------------
    %Create non-georeferenced grid in meters (origin = 0,0)
    x_u   = (x_rho(:,1:L)   + x_rho(:,2:Lp))/2;
    y_u   = (y_rho(:,1:L)   + y_rho(:,2:Lp))/2;

    %Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_u = x_u * cos(rotangle) - y_u * sin(rotangle);
    Ry_u = x_u * sin(rotangle) + y_u * cos(rotangle);

    %Estimate Latitude and Longitude of each Grid point
    lat_u(1,1) = lat + (Ry_u(1,1) ./ latdist);
    lon_u(1,1) = lon + (Rx_u(1,1) ./ spheriq_dist(lon,lat_u(1,1),lon+1,lat_u(1,1)));
    for i = 2:length(Ry_u(:,1));
        lat_u(i,1) = lat + (Ry_u(i,1) ./ latdist);
        lon_u(i,1) = lon_u(i-1,1) + ((Rx_u(i,1)-Rx_u(i-1,1)) ./ spheriq_dist(lon_u(i-1,1),lat_u(i,1),lon_u(i-1,1)+1,lat_u(i,1)));
    end

    for i = 1:length(Ry_u(:,1));
        for j = 2:length(Ry_u(1,:));
            lat_u(i,j) = lat + (Ry_u(i,j) ./ latdist);
            lon_u(i,j) = lon_u(i,j-1) + ((Rx_u(i,j)-Rx_u(i,j-1)) ./ spheriq_dist(lon_u(i,j-1),lat_u(i,j),lon_u(i,j-1)+1,lat_u(i,j)));
        end
    end
    clear x_u y_u Rx_u Ry_u

    % V GRID --------------------------------------------------------------
    %Create non-georeferenced grid in meters (origin = 0,0)
    x_v   = (x_rho(1:M,:)   + x_rho(2:Mp,:))/2;
    y_v   = (y_rho(1:M,:)   + y_rho(2:Mp,:))/2;

    %Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_v = x_v * cos(rotangle) - y_v * sin(rotangle);
    Ry_v = x_v * sin(rotangle) + y_v * cos(rotangle);

    %Estimate Latitude and Longitude of each Grid point
    lat_v(1,1) = lat + (Ry_v(1,1) ./ latdist);
    lon_v(1,1) = lon + (Rx_v(1,1) ./ spheriq_dist(lon,lat_v(1,1),lon+1,lat_v(1,1)));
    for i = 2:length(Ry_v(:,1));
        lat_v(i,1) = lat + (Ry_v(i,1) ./ latdist);
        lon_v(i,1) = lon_v(i-1,1) + ((Rx_v(i,1)-Rx_v(i-1,1)) ./ spheriq_dist(lon_v(i-1,1),lat_v(i,1),lon_v(i-1,1)+1,lat_v(i,1)));
    end

    for i = 1:length(Ry_v(:,1));
        for j = 2:length(Ry_v(1,:));
            lat_v(i,j) = lat + (Ry_v(i,j) ./ latdist);
            lon_v(i,j) = lon_v(i,j-1) + ((Rx_v(i,j)-Rx_v(i,j-1)) ./ spheriq_dist(lon_v(i,j-1),lat_v(i,j),lon_v(i,j-1)+1,lat_v(i,j)));
        end
    end
    clear x_v y_v Rx_v Ry_v

    % PSI GRID ------------------------------------------------------------
    %Create non-georeferenced grid in meters (origin = 0,0)
    x_psi = (x_rho(1:M,1:L) + x_rho(2:Mp,2:Lp))/2;
    y_psi = (y_rho(1:M,1:L) + y_rho(2:Mp,2:Lp))/2;

    %Rotate grid ...See: http://en.wikipedia.org/wiki/Rotation_(mathematics)
    Rx_psi = x_psi * cos(rotangle) - y_psi * sin(rotangle);
    Ry_psi = x_psi * sin(rotangle) + y_psi * cos(rotangle);

    %Estimate Latitude and Longitude of each Grid point
    lat_psi(1,1) = lat + (Ry_psi(1,1) ./ latdist);
    lon_psi(1,1) = lon + (Rx_psi(1,1) ./ spheriq_dist(lon,lat_psi(1,1),lon+1,lat_psi(1,1)));
    for i = 2:length(Ry_psi(:,1));
        lat_psi(i,1) = lat + (Ry_psi(i,1) / latdist);
        lon_psi(i,1) = lon_psi(i-1,1) + ((Rx_psi(i,1)-Rx_psi(i-1,1)) ./ spheriq_dist(lon_psi(i-1,1),lat_psi(i,1),lon_psi(i-1,1)+1,lat_psi(i,1)));
    end

    for i = 1:length(Ry_psi(:,1));
        for j = 2:length(Ry_psi(1,:));
            lat_psi(i,j) = lat + (Ry_psi(i,j) ./ latdist);
            lon_psi(i,j) = lon_psi(i,j-1) + ((Rx_psi(i,j)-Rx_psi(i,j-1)) ./ spheriq_dist(lon_psi(i,j-1),lat_psi(i,j),lon_psi(i,j-1)+1,lat_psi(i,j)));
        end
    end
    clear x_rho y_rho x_psi y_psi Rx_psi Ry_psi i j latdist normradius


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Grid spacing and other grid parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    el      = lat_u(end,1) - lat_u(1,1);
    xl      = lon_v(1,end) - lon_v(1,1);

    %Thanks to John Wilkin for this following section.
    dx=zeros(Mp,Lp);
    dy=zeros(Mp,Lp);

    dx(:,2:L)=spheriq_dist(lon_u(:,2:L),lat_u(:,2:L),lon_u(:,1:Lm),lat_u(:,1:Lm)); %sperical distance calculation
    dx(:,1)=dx(:,2);
    dx(:,Lp)=dx(:,L);
    dy(2:M,:)=spheriq_dist(lon_v(2:M,:),lat_v(2:M,:),lon_v(1:Mm,:),lat_v(1:Mm,:)); %sperical distance calculation
    dy(1,:)=dy(2,:);
    dy(Mp,:)=dy(M,:);
    pm=1./dx;
    pn=1./dy;

    dndx = zeros(size(pm));
    dmde = dndx;
    dndx(2:M,2:L)=0.5*(1./pn(2:M,3:Lp) - 1./pn(2:M,1:Lm));
    dmde(2:M,2:L)=0.5*(1./pm(3:Mp,2:L) - 1./pm(1:Mm,2:L));


    % Coriolis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = 2 .* 7.29E-5 .* sin(lat_rho .* (pi/180)); %Estimation of Coriolis over the grid domain. OMEGA=7.29E-5
    %More info: http://en.wikipedia.org/wiki/Coriolis_effect#Formula


    % Bathymetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isnan(hh);
        h = hh*ones(size(x));
        h = ones(length(y),1)*h(:)';
    else
        load(Bathy);
        h = griddata(xbathy,ybathy,zbathy,lon_rho,lat_rho,'linear');
        h(isnan(h)) = -1;
    end
    h(h<0) = 0;   % Flatten hills and mountains (i.e. positive depths)
    h = h + minh; % Add the depth offset minh (specified in USER SETTINGS)
    clear xbathy ybathy zbathy
    %NOTE: Bathymetry smoothing occurs below "generating masks"




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERATING MASKS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['DONE!',char(13),char(13),char(13),'GENERATING masks...'])
    
    %  Land/Sea mask on RHO-points... and a NaN version of the mask for plotting.
    mask_rho = ones(size(lon_rho));
    mask_rho_nan = mask_rho;

    mask_rho(h <= minh) = 0;
    mask_rho_nan(h <= minh) = nan;

    %  Land/Sea mask on U-points.
    mask_u = mask_rho(:,1:L) .* mask_rho(:,2:Lp);

    %  Land/Sea mask on V-points.
    mask_v = mask_rho(1:M,:) .* mask_rho(2:Mp,:);

    %  Land/Sea mask on PSI-points.
    mask_psi = mask_u(1:M,:) .* mask_u(2:Mp,:);
    disp('DONE!')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp([char(13),char(13)])
    disp('---------------------------------------')
    disp('WORKING with previously-generated GRID!')
    disp('---------------------------------------')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END of grid generation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Check for a grid variable. If it is not there... Warn user to TURN ON create_grid in USER SETTINGS
if exist('lat_rho','var') == 0;
    disp([char(13),char(13)])
    warning('ABSENT GRID VARIABLE!!!')
    disp('NOTE: You may have to SWITCH ON create_grid in the USER SETTINGS section and run EASYGRID again')
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHING Bathymetry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if smooth_grid == 1; % Only smooth bathymetry if switch (in USER SETTINGS) is turned ON
    disp([char(13),char(13),'SMOOTHING Bathymetry...'])
    h = smth_bath(h,mask_rho); %Smooth the grid
    disp('DONE!')
    clear smoothing
else
    clear smoothing
end
%--------------------------------------------------------------------------



  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCREEN DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if screen_output == 1; % Only display on screen if switch (in USER SETTINGS) is turned ON
    disp([char(13),char(13),char(13)])
    disp(['SCREEN DISPLAY:',char(13)])
    disp('COPY-PASTE the following parameters into your ocean.in file')
    disp('----------------------------------------------------------------------------------------------')
    disp(char(13))
    disp(['    Lm == ' num2str(Lm) '         ! Number of I-direction INTERIOR RHO-points'])
    disp(['    Mm == ' num2str(Mm) '         ! Number of J-direction INTERIOR RHO-points'])
    disp(['     N == ' num2str(N) '          ! Number of vertical levels'])
    disp(char(13))
    disp(['Make sure the Baroclinic time-step (DT) in your ocean.in file is less than: ' num2str(sqrt(((min(min(dx))^2)+(min(min(dy))^2)) / (9.8 * (min(min(h))^2)))) ' seconds'])
    disp('----------------------------------------------------------------------------------------------')
end
%--------------------------------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_grid == 1; % Only plot if switch (in USER SETTINGS) is turned ON
    disp([char(13),'PLOTTING grid...'])
    % check if you have John Wilkin's pcolor... if not use the normal pcolor
    cla % Erase axis but keep figure (this is to refresh old plot with new plot, while keeping figure size).
    if exist('pcolorjw.m') == 2;
        pcolorjw(lon_rho,lat_rho,h.*mask_rho_nan);shading faceted;cb=colorbar; title(cb,[{'Depth'};{'(m)'}])
    else
        pcolor(lon_rho,lat_rho,h.*mask_rho_nan);shading faceted;cb=colorbar; title(cb,[{'Depth'};{'(m)'}])
    end
    axis square
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)');
    title(['ROMS grid: ' Grid_filename]);    
    % Check if there is a coastline file
    if exist('Coast') == 1;
        hold on
        plot(Coast.lon,Coast.lat,'k');
    end
    clear cb
    disp(['DONE!'])
end
%--------------------------------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASK EDITING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if edit_mask == 1; % Only edit mask if switch (in USER SETTINGS) is turned ON
    disp([char(13),'EDITING MASK...'])
    disp('To FINISH editing... press right-button (on mouse).')
    figure(2)
    pcolorjw(lon_rho,lat_rho,mask_rho);caxis([0 1]);shading faceted;colormap([0.6 0.4 0;0.6 0.9 1]);colorbar;
    axis square
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)');
    title('LEFT-button: sea/land toggle ... RIGHT-button: Finish!','BackgroundColor','red');
    hold on
    if exist('Coast') == 1; 
        plot(Coast.lon,Coast.lat,'k');
    end
    button = 1;
    while button == 1
        [xi,yi,button] = ginput(1);
        if button == 1;
            costfunct = ((lon_rho-xi).^2) + ((lat_rho-yi).^2);
            [loni,lonj] = find(costfunct == min(min(costfunct)));
                if mask_rho(loni,lonj) == 0;
                    mask_rho(loni,lonj) = 1;
                else
                    mask_rho(loni,lonj) = 0;
                end
            pcolorjw(lon_rho,lat_rho,mask_rho);caxis([0 1]);shading faceted;colormap([0.6 0.4 0;0.6 0.9 1]);colorbar;
                if exist('Coast') == 1;
                    plot(Coast.lon,Coast.lat,'k');
                end
        end
        % Update u, v, psi and rho_nan MASKS
        mask_rho_nan = mask_rho;
        mask_rho_nan(mask_rho == 0) = nan;
        mask_u = mask_rho(:,1:L) .* mask_rho(:,2:Lp);
        mask_v = mask_rho(1:M,:) .* mask_rho(2:Mp,:);
        mask_psi = mask_u(1:M,:) .* mask_u(2:Mp,:); 
    end
    clear button xi yi costfunct
    disp(['DONE!'])
end
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE GRID FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_grid == 1; % Only save grid if switch (in USER SETTINGS) is turned ON
    try
        % creating grid files
        grd_name   = [Grid_filename '_grd.nc'];
        nc         = netcdf(grd_name,'clobber');
        nc.description = Descrip_grd;
        nc.author  = Author;
        nc.created = [ 'by ' nc.author ' on ' Computer 'at '  datestr(now) ', using roms_make_simple_grid.m'];
        nc.type    = 'grid file';
        disp(char(13))
        disp(['SAVING Grid NetCDF file: ' grd_name '   ...'])
        %Note: Other NetCDF definitions are specifies at the beginning of this code under the USER SETTINGS section.

        % Dimensions
        nc('xi_rho') = Lp;
        nc('xi_u')   = L;
        nc('xi_v')   = Lp;
        nc('xi_psi') = L;

        nc('eta_rho')= Mp;
        nc('eta_u')  = Mp;
        nc('eta_v')  = M;
        nc('eta_psi')= M;

        nc('one')    = 1;


        % Create variables to NetCDF
        dims                  = { 'eta_rho'; 'xi_rho'};
        nc{ 'dmde'}           = ncdouble(dims);
        nc{ 'dmde'}(:,:)      = dmde;
        nc{ 'dmde'}.long_name = 'eta derivative of inverse metric factor pm' ;
        nc{ 'dmde'}.units     = 'meter' ;
        %--------------------------------------------------------------------------
        dims                  = { 'eta_rho'; 'xi_rho'};
        nc{ 'dndx'}           = ncdouble(dims);
        nc{ 'dndx'}(:,:)      = dndx;
        nc{ 'dndx'}.long_name = 'xi derivative of inverse metric factor pn' ;
        nc{ 'dndx'}.units     = 'meter' ;
        %--------------------------------------------------------------------------
        dims                  = { 'one'};
        nc{ 'el'}             = ncdouble(dims);
        nc{ 'el'}(:)          = el;
        nc{ 'el'}.long_name   = 'domain length in the ETA-direction' ;
        nc{ 'el'}.units       = 'degrees' ;
        %--------------------------------------------------------------------------
        dims                  = { 'eta_rho'; 'xi_rho'};
        nc{ 'f'}              = ncdouble(dims);
        nc{ 'f'}(:,:)         = f;
        nc{ 'f'}.long_name    = 'Coriolis parameter at RHO-points' ;
        nc{ 'f'}.units        = 'second-1' ;
        %--------------------------------------------------------------------------
        dims                  = { 'eta_rho'; 'xi_rho'};
        nc{ 'h'}              = ncdouble(dims);
        nc{ 'h'}(:,:)         = h;
        nc{ 'h'}.long_name    = 'bathymetry at RHO-points' ;
        nc{ 'h'}.units        = 'meter' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_rho'; 'xi_rho'};
        nc{ 'lat_rho'}           = ncdouble(dims);
        nc{ 'lat_rho'}(:,:)      = lat_rho;
        nc{ 'lat_rho'}.long_name = 'latitude of RHO-points' ;
        nc{ 'lat_rho'}.units     = 'degree_north' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_psi'; 'xi_psi'};
        nc{ 'lat_psi'}           = ncdouble(dims);
        nc{ 'lat_psi'}(:,:)      = lat_psi;
        nc{ 'lat_psi'}.long_name = 'latitude of PSI-points' ;
        nc{ 'lat_psi'}.units     = 'degree_north' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_u'; 'xi_u'};
        nc{ 'lat_u'}             = ncdouble(dims);
        nc{ 'lat_u'}(:,:)        = lat_u;
        nc{ 'lat_u'}.long_name   = 'latitude of U-points' ;
        nc{ 'lat_u'}.units       = 'degree_north' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_v'; 'xi_v'};
        nc{ 'lat_v'}             = ncdouble(dims);
        nc{ 'lat_v'}(:,:)        = lat_v;
        nc{ 'lat_v'}.long_name   = 'latitude of V-points' ;
        nc{ 'lat_v'}.units       = 'degree_north' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_rho'; 'xi_rho'};
        nc{ 'lon_rho'}           = ncdouble(dims);
        nc{ 'lon_rho'}(:,:)      = lon_rho;
        nc{ 'lon_rho'}.long_name = 'longitude of RHO-points' ;
        nc{ 'lon_rho'}.units     = 'degree_east' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_psi'; 'xi_psi'};
        nc{ 'lon_psi'}           = ncdouble(dims);
        nc{ 'lon_psi'}(:,:)      = lon_psi;
        nc{ 'lon_psi'}.long_name = 'longitude of PSI-points' ;
        nc{ 'lon_psi'}.units     = 'degree_east' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_u'; 'xi_u'};
        nc{ 'lon_u'}             = ncdouble(dims);
        nc{ 'lon_u'}(:,:)        = lon_u;
        nc{ 'lon_u'}.long_name   = 'longitude of U-points' ;
        nc{ 'lon_u'}.units       = 'degree_east' ;
        %--------------------------------------------------------------------------
        dims                     = { 'eta_v'; 'xi_v'};
        nc{ 'lon_v'}             = ncdouble(dims);
        nc{ 'lon_v'}(:,:)        = lon_v;
        nc{ 'lon_v'}.long_name   = 'longitude of V-points' ;
        nc{ 'lon_v'}.units       = 'degree_east' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_rho'; 'xi_rho'};
        nc{ 'mask_rho'}           = ncdouble(dims);
        nc{ 'mask_rho'}(:,:)      = mask_rho;
        nc{ 'mask_rho'}.long_name = 'mask on RHO-points' ;
        nc{ 'mask_rho'}.option_0  = 'land' ;
        nc{ 'mask_rho'}.option_1  = 'water' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_psi'; 'xi_psi'};
        nc{ 'mask_psi'}           = ncdouble(dims);
        nc{ 'mask_psi'}(:,:)      = mask_psi;
        nc{ 'mask_psi'}.long_name = 'mask on PSI-points' ;
        nc{ 'mask_psi'}.option_0  = 'land' ;
        nc{ 'mask_psi'}.option_1  = 'water' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_u'; 'xi_u'};
        nc{ 'mask_u'}             = ncdouble(dims);
        nc{ 'mask_u'}(:,:)        = mask_u;
        nc{ 'mask_u'}.long_name   = 'mask on U-points' ;
        nc{ 'mask_u'}.option_0    = 'land' ;
        nc{ 'mask_u'}.option_1    = 'water' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_v'; 'xi_v'};
        nc{ 'mask_v'}             = ncdouble(dims);
        nc{ 'mask_v'}(:,:)        = mask_v;
        nc{ 'mask_v'}.long_name   = 'mask on V-points' ;
        nc{ 'mask_v'}.option_0    = 'land' ;
        nc{ 'mask_v'}.option_1    = 'water' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_rho'; 'xi_rho'};
        nc{ 'pm'}                 = ncdouble(dims);
        nc{ 'pm'}(:,:)            = pm;
        nc{ 'pm'}.long_name       = 'curvilinear coordinate metric in XI' ;
        nc{ 'pm'}.units           = 'meter-1' ;
        %--------------------------------------------------------------------------
        dims                      = { 'eta_rho'; 'xi_rho'};
        nc{ 'pn'}                 = ncdouble(dims);
        nc{ 'pn'}(:,:)            = pn;
        nc{ 'pn'}.long_name       = 'curvilinear coordinate metric in ETA' ;
        nc{ 'pn'}.units           = 'meter-1' ;
        %--------------------------------------------------------------------------
        dims                      = { 'one'};
        nc{ 'spherical'}          = ncchar(dims);
        nc{ 'spherical'}(:)       = 'F';
        nc{ 'spherical'}.long_name= 'Grid type logical switch' ;
        nc{ 'spherical'}.option_T = 'spherical' ;
        %--------------------------------------------------------------------------
        dims                      = {'one'};
        nc{ 'xl'}                 = ncdouble(dims);
        nc{ 'xl'}(:)              = xl;
        nc{ 'xl'}.long_name       = 'domain length in the XI-direction' ;
        nc{ 'xl'}.units           = 'degrees' ;
        %--------------------------------------------------------------------------
        %--------------------------------------------------------------------------
        % Parameters
        dims                  = { 'one'};
        nc{ 'X'}              = ncdouble(dims);
        nc{ 'X'}(:)           = X;
        nc{ 'X'}.description  = 'width of domain (degrees)';
        %--------------------------------------------------------------------------
        dims                  = { 'one'};
        nc{ 'Y'}              = ncdouble(dims);
        nc{ 'Y'}(:)           = Y;
        nc{ 'Y'}.description  = 'length of domain (degrees)';
        %--------------------------------------------------------------------------
        dims                  = { 'one'};
        nc{ 'dx'}             = ncdouble(dims);
        nc{ 'dx'}(:)          = (resol ./ (60*1852)) .* spheriq_dist(lon,lat,lon+1,lat);;   % Estimated resolution in x (degrees)0.002;
        nc{ 'dx'}.description = 'resolution in x (degrees)';
        %--------------------------------------------------------------------------
        dims                  = { 'one'};
        nc{ 'dy'}             = ncdouble(dims);
        nc{ 'dy'}(:)          = resol ./ (60*1852);         % Estimated resolution in y (degrees)
        nc{ 'dy'}.description = 'resolution in y (degrees)';
        
        disp('DONE!')
        close(nc);
        clear nc
    catch
        warning('Problem writing .nc grid file. You may not have MEXCDF installed...')
    end
end
%--------------------------------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE INITIAL CONDITIONS FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_init == 1; % Only save init_file if switch (in USER SETTINGS) is turned ON
    try
        % creating initial cond file
        ini_name   = [Grid_filename '_ini.nc'];
        nc         = netcdf(ini_name,'clobber');

        nc.description = Descrip_ini;
        nc.author  = Author;
        nc.created = [ 'by ' nc.author ' on ' Computer 'at '  datestr(now) ', using roms_make_simple_grid.m'];
        nc.type    = 'grid file';
        disp(char(13))
        disp(['SAVING Initial Conditions NetCDF file: ' grd_name '   ...'])
        %Note: Other NetCDF definitions are specifies at the beginning of this
        %code under the USER SETTINGS section.

        % Dimensions
        nc('xi_rho') = Lp;
        nc('xi_u')   = L;
        nc('xi_v')   = Lp;
        nc('xi_psi') = L;

        nc('eta_rho')= Mp;
        nc('eta_u')  = Mp;
        nc('eta_v')  = M;
        nc('eta_psi')= M;

        nc('s_rho')  = N;
        nc('time')   = 1;

        % Create variables to NetCDF
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'NH4'}                 = ncdouble(dims);
        nc{ 'NH4'}(:,:,:)          = ones(N,Mp,Lp).* NH4;
        nc{ 'NH4'}.time            = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'NO3'}                 = ncdouble(dims);
        nc{ 'NO3'}(:,:,:)          = ones(N,Mp,Lp).* NO3;
        nc{ 'NO3'}.time            = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'chlorophyll1'}        = ncdouble(dims);
        nc{ 'chlorophyll1'}(:,:,:) = ones(N,Mp,Lp).* chlorophyll1;
        nc{ 'chlorophyll1'}.time   = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'chlorophyll2'}        = ncdouble(dims);
        nc{ 'chlorophyll2'}(:,:,:) = ones(N,Mp,Lp).* chlorophyll2;
        nc{ 'chlorophyll2'}.time   = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'detritus1'}           = ncdouble(dims);
        nc{ 'detritus1'}(:,:,:)    = ones(N,Mp,Lp).* detritus1;
        nc{ 'detritus1'}.time      = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'detritus2'}           = ncdouble(dims);
        nc{ 'detritus2'}(:,:,:)    = ones(N,Mp,Lp).* detritus2;
        nc{ 'detritus2'}.time      = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'detritusC1'}          = ncdouble(dims);
        nc{ 'detritusC1'}(:,:,:)   = ones(N,Mp,Lp).* detritusC1;
        nc{ 'detritusC1'}.time     = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'detritusC2'}          = ncdouble(dims);
        nc{ 'detritusC2'}(:,:,:)   = ones(N,Mp,Lp).* detritusC2;
        nc{ 'detritusC2'}.time     = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'};
        nc{ 'ocean_time'}          = ncdouble(dims);
        nc{ 'ocean_time'}(:)       = 53004;
        nc{ 'ocean_time'}.units    = 'days since 1858-11-17 00:00:00 GMT' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'phyto1'}              = ncdouble(dims);
        nc{ 'phyto1'}(:,:,:)       = ones(N,Mp,Lp).* phyto1;
        nc{ 'phyto1'}.time         = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'phyto2'}              = ncdouble(dims);
        nc{ 'phyto2'}(:,:,:)       = ones(N,Mp,Lp).* phyto2;
        nc{ 'phyto2'}.time         = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'phytoC1'}             = ncdouble(dims);
        nc{ 'phytoC1'}(:,:,:)      = ones(N,Mp,Lp).* phytoC1;
        nc{ 'phytoC1'}.time        = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'phytoC2'}             = ncdouble(dims);
        nc{ 'phytoC2'}(:,:,:)      = ones(N,Mp,Lp).* phytoC2;
        nc{ 'phytoC2'}.time        = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'salt'}                = ncdouble(dims);
        nc{ 'salt'}(:,:,:)         = ones(N,Mp,Lp).* salt;
        nc{ 'salt'}.time           = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'temp'}                = ncdouble(dims);
        nc{ 'temp'}(:,:,:)         = ones(N,Mp,Lp).* temp;
        nc{ 'temp'}.time           = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_u'; 'xi_u'};
        nc{ 'u'}                   = ncdouble(dims);
        nc{ 'u'}(:,:,:)            = ones(N,Mp,L).* u;
        nc{ 'u'}.time              = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 'eta_u'; 'xi_u'};
        nc{ 'ubar'}                = ncdouble(dims);
        nc{ 'ubar'}(:,:,:)         = ones(Mp,L).* ubar;
        nc{ 'ubar'}.time           = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_v'; 'xi_v'};
        nc{ 'v'}                   = ncdouble(dims);
        nc{ 'v'}(:,:,:)            = ones(N,M,Lp).* v;
        nc{ 'v'}.time              = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 'eta_v'; 'xi_v'};
        nc{ 'vbar'}                = ncdouble(dims);
        nc{ 'vbar'}(:,:,:)         = ones(M,Lp).* vbar;
        nc{ 'vbar'}.time           = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 'eta_rho'; 'xi_rho'};
        nc{ 'zeta'}                = ncdouble(dims);
        nc{ 'zeta'}(:,:,:)         = ones(Mp,Lp).* zeta;
        nc{ 'zeta'}.time           = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'zooplankton'}        = ncdouble(dims);
        nc{ 'zooplankton'}(:,:,:) = ones(N,Mp,Lp).* zooplankton;
        nc{ 'zooplankton'}.time   = 'ocean_time' ;
        %--------------------------------------------------------------------------
        dims                       = { 'time'; 's_rho'; 'eta_rho'; 'xi_rho'};
        nc{ 'zooplanktonC'}        = ncdouble(dims);
        nc{ 'zooplanktonC'}(:,:,:) = ones(N,Mp,Lp).* zooplanktonC;
        nc{ 'zooplanktonC'}.time   = 'ocean_time' ;
        %--------------------------------------------------------------------------

        disp('DONE!')
        close(nc);

    catch
        warning('Problem writing .nc initial-conditions file. You may not have MEXCDF installed...')
    end
end
disp(char(13))
disp('***************************************************************')
disp(['elapsed time: ' num2str(toc/60) ' minute(s)'])
disp('***************************************************************')
disp('EASYGRID is DONE!!! *******************************************')
disp('***************************************************************')
disp([char(13),char(13),char(13)])

clear nc Author Bathy Computer Descrip_grd Descrip_ini Grid_filename X Y dims ...
    grd_name ini_name plot_grid polar_rad npass order rlim ...
    create_grid edit_mask screen_output save_grid save_init smooth_grid
