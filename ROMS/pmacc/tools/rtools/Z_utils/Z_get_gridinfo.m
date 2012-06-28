function ginfo = Z_get_gridinfo(basename,bathy_dir)
% OBSOLETE 10/21/2011 PM
%-------------------------------------------------------------
% Z_get_gridinfo.m  6/27/2007  Parker MacCready
%
% This is where you:
%
% (i) define all the characteristics of a grid, and
%
% (ii) return them to the calling function in a structure "ginfo"
%
% NO ARGUMENTS: list all possibilites to the screen
% ARGUMENTS = (#,0) return limited "ginfo" for basename number #
% ARGUMENTS = (#,1) return full "ginfo" for basename number #
%
% bathy_dir is the directory with the bathymetry in it
%
% edited by DAS, 2009 to use SNCTOOLS, added grids, salish for example uses
% a different stretching technique (parabola)
% SNG added ptx_high Mar 2011 and also included enhanced plotting options
% also changed input to include bdir such that paths can all be set in
% make_grid to reduce changes when one switches from a local to a remote
% machine
%-------------------------------------------------------------


info_switch = 1;
if nargin == 0
    info_switch = 0;
    bn_num_list = 1:100; % OK to use tenths
    disp(' ')
else
    bn_num_list = 1;
    %bn_num_list = runnumber;
end

for jj = 1:length(bn_num_list)
    bn_num = bn_num_list(jj);
    disp_on = logical(1);
    switch basename
    	case 'ahab7a'
    	    description = 'grid for PS-AHAB';
        	if info_switch==1
			    % specify grid resolution at some anchor points
				lonmin = -129;
				lonmax = -122.05;
				latmin = 43;
				latmax = 51.5;
				x =       [lonmin  -125  -123.1  -122.3  lonmax];
				resx_km = [     6     3     0.8     0.8     1.6];
				y =       [latmin   45.8  46.25  46.5  47.2  48.1   49  latmax];
				resy_km = [     6      3    0.8     3   1.2   1.2  1.6       6];
				r = 1.1;				
				
                km2m = 1000; %convert km to m
                RE = Z_RE(latmin,latmax); %earth radius at mean lat
                latm = (latmin+latmax)/2; %mean latitude
                %dl = 1/80;
                hmin = 4;
                sph = 'T';
			
        		% make the lon-lat grid
        		resx = resx_km ./ (pi*RE/km2m/180) ./ cosd(latm);
        		resy = resy_km ./ (pi*RE/km2m/180);
        		lonr = resolution2grid(x,resx,r);
        		latr = resolution2grid(y,resy,r);
        		
        		% load bathymetry and interpolate
                bname = [bathy_dir,'hybrid/pnw_combined_full.mat'];
                load(bname);
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
                
                % estimate speed relative to ptx_high
                demand = prod(size(Lonr)) / (570*249) / min( min(resy_km), min(resx_km) );
    			disp(['run time = ' num2str(demand) ' relative to ptx_high']);
    		end
    	
        case 'ptx_high'
            description = 'stretched, high-res grid for PNWTOX';
            if info_switch == 1 % only do this if we are creating this grid
                % higher res ps grid, try swath approach instead of parabola
                %set grid limits
                lonmin = -127;
                lonmax = -122;
                latmin = 43;
                latmax = 50;
                km2m = 1000; %convert km to m
                RE = Z_RE(latmin,latmax); %earth radius at mean lat
                latm = (latmin+latmax)/2; %mean latitude
                %dl = 1/80;
                hmin = 4;
                sph = 'T';
                %set min and max resolutions (km) in each tanh function
                minresw = 1; 
                maxresw = 3;
                minrese = minresw; 
                maxrese = 2;
                minress = 1;
                maxress = 2;
                minresn = minress;
                maxresn = 2;
                %to set up the tanh profiles, need to calculate delta_res
                dresw = maxresw-minresw;
                drese = maxrese-minrese;
                dress = maxress-minress;
                dresn = maxresn-minresn;
                %set scale factors to determine the width over which the
                %tanh decays (higher number = wider tanh)
                %widens it to about 2*scalew km on either side, for example
                %scalew = 40 makes a tanh nearly 200km wide, scalew = 20
                %makes a tanh nearly 100km wide
                scalew = 40; 
                scalee = 30;
                scales = 50;
                scalen = 40;
                % x is km with zero at lower left corner
                % y is km with zero at lower left corner
                lonr_0 = linspace(lonmin,lonmax,200);
                latr_0 = linspace(latmin,latmax,200);
                x=(pi*RE*cosd(latm)/180)*(lonr_0-(lonmin))/km2m;
                y=(pi*RE/180)*(latr_0-(latmin))/km2m;
                xt1 = x(dsearchn(lonr_0',-125.7)); %target lon to start tanh decay at
                yt1 = y(dsearchn(latr_0',45.5));
                xt2 = x(dsearchn(lonr_0',-123.25)); %target lon to end tanh decay at
                yt2 = y(dsearchn(latr_0',49));
                % shape is a function to multiply dl by
                shape_w = ((minresw+dresw/2)+(dresw/2).*tanh(-(x-xt1)/scalew));
                shape_e = ((minrese+drese/2)+(drese/2).*tanh((x-xt2)/scalee));
                shape = max(shape_e,shape_w);
                %shape = shape_w;
                %
                % DEBUGGING
                plot_shape = 1;
                if plot_shape
                    figure
                    subplot(211)
                    plot(x,shape)
                    title('EW Shape')
                    grid on
                    xlim([0 round(max(y))])
                    ylabel('resolution (km)')
                end
                i = 1;
                lonr(i) = lonmin;
                while lonr(i)<=lonmax
                    i = i+1;
                    shape_fact = interp1(lonr_0,shape,lonr(i-1));
                    lonr(i)=lonr(i-1)+km2m*shape_fact*(180/(pi*RE*cosd(latm)));
                    %lonr(i)=lonr(i-1)+dl*shape_fact;
                end
                title(['EW Shape, scale = ' num2str(scalew) ', ' num2str(scalee) ' N_x = ' num2str(length(lonr))])
                
                % Get the latitude for an isotropic grid
                shape_s = ((minress+dress/2)+(dress/2).*tanh(-(y-yt1)/scales));
                shape_n = ((minresn+dresn/2)+(dresn/2).*tanh((y-yt2)/scalen));
                shape = max(shape_n,shape_s);
                % increase resolution a bit in the CR channel
                %shape = shape .* (1 - .6*exp(-((y-290).*(y-290))/(30*30)));

                % DEBUGGING
                if plot_shape
                    subplot(212)
                    plot(y,shape)
                    grid on
                    xlim([0 round(max(y))])
                    xlabel('distance (km)')
                end
                set(gcf,'position',[350 1000 550 400]);

                i=1;
                latr(i)=latmin;
                while latr(i)<=latmax
                    i=i+1;
                    shape_fact = interp1(latr_0,shape,latr(i-1));
                    latr(i) = latr(i-1)+km2m*shape_fact*(180/(pi*RE));
                    %latr(i) = latr(i-1)+km2m*shape_fact*(180/(pi*RE))*cosd(latr(i-1));
                end
                title(['NS Shape, scale = ' num2str(scalen) ', ' num2str(scales) ' N_y = ' num2str(length(latr))])

                % load bathymetry and interpolate
                bname = [bathy_dir,'hybrid/pnw_combined_full.mat'];
                load(bname);
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
            end
        case 'salish'
            description = 'stretched PS grid';
            if info_switch == 1 % only do this if we are creating this grid
                %
                % like first davetest grid, but with finer resolution in PS
                % and extending the fine reslution part north/south and
                % east/west, using parabola like fxn to stretch grid
                % instead of tanh
                %
                lonmin = -127;
                lonmax = -122;
                latmin = 45;
                latmax = 50;
                dl = 1/270;
                upscale = 11; %multiply dl by in the shape fxn below
                hmin = 4;
                sph = 'T';
                % x is km with zero at lower left corner
                % y is km with zero at lower left corner
                mlr = pi*((latmin+latmax)/2)/180;
                lonr_0 = linspace(lonmin,lonmax,200);
                latr_0 = linspace(latmin,latmax,200);
                x=111.3*cos(mlr)*(lonr_0-(lonmin));
                y=111.3*(latr_0-(latmin));
                xt = x(dsearchn(lonr_0',-122.7)); %target lon for center of parabola
                xt1 = dsearchn(lonr_0',-125.8); %target index to end parabola at to west
                yt = y(dsearchn(latr_0',47.7)); %target lat for center of parabola
                yt1 = dsearchn(latr_0',45.5); %target index to end parabola at to S
                yt2 = dsearchn(latr_0',50.5); %target index to end parabola at to N
                % shape is a function to multiply dl by   
                shape = upscale*ones(size(x));
                shape_p = 1+(upscale-1)*((x(xt1:end)-xt).^2)./max((x(xt1:end)-xt).^2);
                shape(xt1:end) = shape_p;
                % DEBUGGING
                plot_shape = 1;
                if plot_shape
                    figure
                    subplot(211)
                    plot(x,shape)
                    title('EW Shape')
                    grid on
                end
                i = 1;
                lonr(i) = lonmin;
                while lonr(i)<=lonmax
                    i = i+1;
                    shape_fact = interp1(lonr_0,shape,lonr(i-1));
                    lonr(i)=lonr(i-1)+dl*shape_fact;
                end
                % Get the latitude for an isotropic grid
                shape = upscale*ones(size(y));
                shape_p = 1+(upscale-1)*((y(yt1:yt2)-yt).^2)./max((y(yt1:yt2)-yt).^2);
                shape(yt1:yt2) = shape_p;
                % increase resolution a bit in the CR channel
                shape = shape .* (1 - .6*exp(-((y-135).*(y-135))/(15*15)));
                % DEBUGGING
                if plot_shape
                    subplot(212)
                    plot(y,shape)
                    title('NS Shape')
                    grid on
                end
                i=1;
                latr(i)=latmin;
                while latr(i)<=latmax
                    i=i+1;
                    shape_fact = interp1(latr_0,shape,latr(i-1));
                    latr(i)=latr(i-1)+dl*shape_fact*cos(latr(i-1)*pi/180);
                end
                % load bathymetry and interpolate, this uses the Cascadia
                % bathymetry dataset
                %
                %bdir='/Users/daves/Desktop/PugetSound/R3_tools/In_data/Topo/cascadia_topo/';
                %bname = [bdir,'cascadia_gridded.mat'];
                
                bdir='/Users/sarahgid/Documents/Research/RTOOLS/In_data/Topo/';
                bname = [bdir,'hybrid/pnw_combined_full.mat'];
                                
                load(bname);           
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
                if(plot_shape)
                    figure;plot(Lonr,Latr,'r.')
                    hold on; %plot_WAcoast;
                    axis image; axis([lonmin lonmax latmin latmax])
                    latind = dsearchn(Latr(:,1),48.3); %line through sjdf
                    dx=sw_dist(Latr(latind,:),Lonr(latind,:),'km');
                    figure; subplot(121);
                    plot(Lonr(latind,2:end),dx,'b');hold on;grid on;
                    latind = dsearchn(Latr(:,1),46.2); %line near cr
                    dx=sw_dist(Latr(latind,:),Lonr(latind,:),'km');
                    plot(Lonr(latind,2:end),dx,'r');hold on;
                    legend('48.3 N','46.2 N','location','southwest');
                    subplot(122)
                    lonind = dsearchn(Lonr(1,:)',-122.45); %line down ps
                    dy=sw_dist(Latr(:,lonind),Lonr(:,lonind),'km');
                    plot(Latr(2:end,lonind),dy,'b');hold on;
                    lonind = dsearchn(Lonr(1,:)',-124.7); %line down ps
                    dy=sw_dist(Latr(:,lonind),Lonr(:,lonind),'km');
                    plot(Latr(2:end,lonind),dy,'r');hold on;
                    legend('122.45 W','124.7 W','location','southwest');
                end
            end
    case 'psmid'
            description = 'mid-res PS grid';
            if info_switch == 1 % only do this if we are creating this grid
                % higher res ps grid than pslow, try swath approach instead of parabola 
                lonmin = -126;
                lonmax = -122;
                latmin = 45;
                latmax = 50;
                dl = 1/150;
                hmin = 4;
                sph = 'T';
                upscale = 8;
                % x is km with zero at lower left corner
                % y is km with zero at lower left corner
                mlr = pi*((latmin+latmax)/2)/180;
                lonr_0 = linspace(lonmin,lonmax,200);
                latr_0 = linspace(latmin,latmax,200);
                x=111.3*cos(mlr)*(lonr_0-(lonmin));
                y=111.3*(latr_0-(latmin));
                xt1 = x(dsearchn(lonr_0',-123.7)); %target lon to start tanh decay at
                yt1 = y(dsearchn(latr_0',46.9));
                xt2 = x(dsearchn(lonr_0',-122)); %target lon to end tanh decay at
                yt2 = y(dsearchn(latr_0',48.75));
                % shape is a function to multiply dl by
                shape_w = (1+upscale*(1+tanh(-(x-xt1)/30))/2);
                shape_e = (1+upscale*(1+tanh((x-xt2)/15))/2);
                %shape = max(shape_e,shape_w);
                shape = shape_w;
                % increase resolution a bit in the CR channel
                %shape = shape .* (1 - .4*exp(-((x-7).*(x-7))/(20*20)));
                % DEBUGGING
                plot_shape = 1;
                if plot_shape
                    figure
                    subplot(211)
                    plot(x,shape)
                    title('EW Shape')
                    grid on
                end
                i = 1;
                lonr(i) = lonmin;
                while lonr(i)<=lonmax
                    i = i+1;
                    shape_fact = interp1(lonr_0,shape,lonr(i-1));
                    lonr(i)=lonr(i-1)+dl*shape_fact;
                end
                % Get the latitude for an isotropic grid
                shape_n = (1+upscale*(1+tanh(-(y-yt1)/15))/2);
                shape_s = (1+upscale*(1+tanh((y-yt2)/30))/2);
                shape = max(shape_n,shape_s);
                %shape = shape_n;
                % DEBUGGING
                if plot_shape
                    subplot(212)
                    plot(y,shape)
                    title('NS Shape')
                    grid on
                end
                i=1;
                latr(i)=latmin;
                while latr(i)<=latmax
                    i=i+1;
                    shape_fact = interp1(latr_0,shape,latr(i-1));
                    latr(i)=latr(i-1)+dl*shape_fact*cos(latr(i-1)*pi/180);
                end
                % load bathymetry and interpolate
                bdir='/Users/daves/Desktop/PugetSound/R3_tools/In_data/Topo/cascadia_topo/';
                bname = [bdir,'cascadia_gridded.mat'];
                load(bname);                        
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
                if(plot_shape)
                    figure;plot(Lonr,Latr,'r.')
                    hold on; plot_WAcoast;
                    axis image
                    k=dsearchn(Latr(:,1),48);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'k');
                    k=dsearchn(Latr(:,1),46.5);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'b');
                end
            end               
   case 'pslow'
            description = 'low-res PS grid';
            if info_switch == 1 % only do this if we are creating this grid
                % low res ps grid, for faster computing and debugging 
                lonmin = -126;
                lonmax = -122;
                latmin = 45;
                latmax = 50;
                dl = 1/60;
                hmin = 3;
                sph = 'T';
                % x is km with zero at lower left corner
                % y is km with zero at lower left corner
                mlr = pi*((latmin+latmax)/2)/180;
                lonr_0 = linspace(lonmin,lonmax,200);
                latr_0 = linspace(latmin,latmax,200);
                x=111.3*cos(mlr)*(lonr_0-(lonmin));
                y=111.3*(latr_0-(latmin));
                xt1 = x(dsearchn(lonr_0',-123.7)); %target lon to start tanh decay at
                yt1 = y(dsearchn(latr_0',46.9));
                xt2 = x(dsearchn(lonr_0',-122)); %target lon to end tanh decay at
                yt2 = y(dsearchn(latr_0',48.75));
                % shape is a function to multiply dl by
                shape_w = (1+7*(1+tanh(-(x-xt1)/30))/2);
                shape_e = (1+7*(1+tanh((x-xt2)/15))/2);
                %shape = max(shape_e,shape_w);
                shape = shape_w;
                % increase resolution a bit in the CR channel
                %shape = shape .* (1 - .4*exp(-((x-7).*(x-7))/(20*20)));
                % DEBUGGING
                plot_shape = 1;
                if plot_shape
                    figure
                    subplot(211)
                    plot(x,shape)
                    title('EW Shape')
                    grid on
                end
                i = 1;
                lonr(i) = lonmin;
                while lonr(i)<=lonmax
                    i = i+1;
                    shape_fact = interp1(lonr_0,shape,lonr(i-1));
                    lonr(i)=lonr(i-1)+dl*shape_fact;
                end
                % Get the latitude for an isotropic grid
                shape_n = (1+12*(1+tanh(-(y-yt1)/15))/2);
                shape_s = (1+7*(1+tanh((y-yt2)/30))/2);
                shape = max(shape_n,shape_s);
                %shape = shape_n;
                % DEBUGGING
                if plot_shape
                    subplot(212)
                    plot(y,shape)
                    title('NS Shape')
                    grid on
                end
                i=1;
                latr(i)=latmin;
                while latr(i)<=latmax
                    i=i+1;
                    shape_fact = interp1(latr_0,shape,latr(i-1));
                    latr(i)=latr(i-1)+dl*shape_fact*cos(latr(i-1)*pi/180);
                end
                % load bathymetry and interpolate
                bdir='/Users/daves/Desktop/PugetSound/R3_tools/In_data/Topo/cascadia_topo/';
                bname = [bdir,'cascadia_gridded.mat'];
                load(bname);              
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
                if(plot_shape)
                    cfile = '~/Desktop/PugetSound/R3_tools/In_data/Topo/coastlines/coastline_regional.mat';
                    load(cfile);
                    figure;plot(Lonr,Latr,'r.')
                    hold on; plot(lon_coast,lat_coast,'-k')
                    axis image
                    k=dsearchn(Latr(:,1),48);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'k');
                    k=dsearchn(Latr(:,1),46.5);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'b');
                end
            end    
      case 'pnwtox1'
            description = 'stretched, mid-res grid for PNWTOX';
            if info_switch == 1 % only do this if we are creating this grid
                % higher res ps grid, try swath approach instead of parabola 
                lonmin = -128;
                lonmax = -122;
                latmin = 43;
                latmax = 50.5;
                dl = 1/80;
                hmin = 4;
                sph = 'T';
                upscale = 5.5;
                % x is km with zero at lower left corner
                % y is km with zero at lower left corner
                mlr = pi*((latmin+latmax)/2)/180;
                lonr_0 = linspace(lonmin,lonmax,200);
                latr_0 = linspace(latmin,latmax,200);
                x=111.3*cos(mlr)*(lonr_0-(lonmin));
                y=111.3*(latr_0-(latmin));
                xt1 = x(dsearchn(lonr_0',-125.9)); %target lon to start tanh decay at
                yt1 = y(dsearchn(latr_0',45.2));
                xt2 = x(dsearchn(lonr_0',-122)); %target lon to end tanh decay at
                yt2 = y(dsearchn(latr_0',49));
                % shape is a function to multiply dl by
                shape_w = (1+upscale*(1+tanh(-(x-xt1)/80))/2);
                shape_e = (1+upscale*(1+tanh((x-xt2)/15))/2);
                %shape = max(shape_e,shape_w);
                shape = shape_w;
                %
                % DEBUGGING
                plot_shape = 1;
                if plot_shape
                    figure
                    subplot(211)
                    plot(x,shape)
                    title('EW Shape')
                    grid on
                end
                i = 1;
                lonr(i) = lonmin;
                while lonr(i)<=lonmax
                    i = i+1;
                    shape_fact = interp1(lonr_0,shape,lonr(i-1));
                    lonr(i)=lonr(i-1)+dl*shape_fact;
                end
                % Get the latitude for an isotropic grid
                shape_s = (1+upscale*(1+tanh(-(y-yt1)/80))/2);
                shape_n = (1+upscale*(1+tanh((y-yt2)/40))/2);
                shape = max(shape_n,shape_s);
                % increase resolution a bit in the CR channel
                %shape = shape .* (1 - .6*exp(-((y-290).*(y-290))/(30*30)));
                %shape = shape_n;
                % DEBUGGING
                if plot_shape
                    subplot(212)
                    plot(y,shape)
                    title('NS Shape')
                    grid on
                end
                i=1;
                latr(i)=latmin;
                while latr(i)<=latmax
                    i=i+1;
                    shape_fact = interp1(latr_0,shape,latr(i-1));
                    latr(i)=latr(i-1)+dl*shape_fact*cos(latr(i-1)*pi/180);
                end
                % load bathymetry and interpolate
                bdir='/Users/daves/Desktop/PugetSound/R3_tools/In_data/Topo/cascadia_topo/';
                bname = [bdir,'cascadia_gridded.mat'];
                load(bname);                        
                [Lonr,Latr] = meshgrid(lonr,latr);
                h = interp2(lon_topo,lat_topo,-z_topo,Lonr,Latr);
                if(plot_shape)
                    figure;plot(Lonr,Latr,'r.')
                    hold on; plot_WAcoast;
                    axis image
                    k=dsearchn(Latr(:,1),48);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'k');
                    k=dsearchn(Latr(:,1),46.5);
                    dx=sw_dist(Latr(k,:),Lonr(k,:),'km');
                    figure; plot(Lonr(k,2:end),dx,'b');
                end
            end               
        otherwise
            disp_on = logical(0);
  end %end switch statement
  
    if disp_on & nargin==0
        disp(['  ',num2str(bn_num),': <',basename,'> ', ...
            description]);
    end
end

% pack everything into the structure "ginfo"
ginfo.basename = basename;
ginfo.description = description;
if info_switch == 1
    ginfo.Lonr = Lonr; ginfo.Latr = Latr;
    ginfo.h = h;
    ginfo.hmin = hmin;
    ginfo.sph = sph;
end
