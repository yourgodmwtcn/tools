function cruiseDate = Z_PS_iniCruise(r_dir, S, obsFile, polyFile, iniDate, ininame, timeBuffer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z_PS_iniCruise.m
%
% m-file for reading in obsfile with CTD data and inserting it into
%  the initial condition already made (T and S only)
%
% uses polygons for regions in puget sound
%
% outputs mean date of cruise, to compare with IC date
%
% written by DAS 3/9/2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% first read in grid info from ocean_ini.nc
gfile = [r_dir, 'grid.nc']; %assumes already made grid.nc in the run folder
outfile = [r_dir,'Ocn/',ininame];

h = nc_varget(gfile,'h'); % get grid data 
[M,L] = size(h);
mask_rho = nc_varget(gfile,'mask_rho');
lat_rho = nc_varget(gfile,'lat_rho');
lon_rho = nc_varget(gfile,'lon_rho');
[z_rho,z_w] = Z_s2z(h,0*h,S);
zz = reshape(z_rho,S.N,M*L); zz = -zz; %makes positive down
% read in initial condition
salt_roms = nc_varget(outfile,'salt'); salt_old = salt_roms; %for debugging plot below
temp_roms = nc_varget(outfile,'temp'); temp_old = temp_roms; %for debugging plot below

%choose year of cruise and file to work on
load(polyFile); %needs to have polys struct defined
polynames = fieldnames(polys); npoly=length(polynames);

%%% now prepare to use obs_extract.m
tspan = [iniDate - timeBuffer iniDate + timeBuffer]; %time buffer around initialization date
vars = [{'salinity'};{'temperature'}];

tic
% for each polygon in puget sound, filter points to include only those in polygon,
%     then use those to overwrite T/S in ocean_ini fields within that same region
cruiseDate = [];
concat = 1;
for n = 1:npoly
    region = char(polynames(n));
    disp([' Working on ', region, ' region'])
    rpos = getfield(polys,region); rx = rpos(:,1); ry = rpos(:,2);
    try 
        obs = obs_extract(obsFile, vars, tspan, [ry rx]);
    catch
        disp(['no data available in this time span for ' region ' region!']); 
        continue
    end
    % now find points in grid and in cruise data inside this region
    inGrid = inpolygon(lon_rho, lat_rho, rx, ry);
    indgrid = find(inGrid==1);
    if(isempty(obs.time))
        disp(['No stations in the ',region,' region!']); 
    else
        [b,k,m]=unique(obs.castid,'first');
        [b,k2,m]=unique(obs.castid);
        ind = [k k2];
        lonstat = obs.longitude(ind(:,1));
        latstat = obs.latitude(ind(:,1));
        [bb,kk,mm] = unique([lonstat latstat],'rows');
        lonstat=lonstat(kk); latstat=latstat(kk);
        ind = ind(kk,:);
        nstations = length(ind(:,1));
        longrid = lon_rho(inGrid);
        latgrid = lat_rho(inGrid);
        if(nstations <= 2) %can't do TRI <= 2 stations
            dd = 1/250; % duplicate station right next to it by dd degrees
            lonstat(end+1:end+2) = [lonstat(end)+dd; lonstat(end)-dd];
            latstat(end+1:end+2) = [latstat(end); latstat(end)-dd];
            ind = [ind; ind(end,:); ind(end,:)];
        end
        TRI = delaunay(lonstat,latstat);
        if concat==1
         lonstats = lonstat(:); latstats = latstat(:); cruiseDate = mean(obs.time);
         sst = obs.temperature(ind(:,1)); sss = obs.salinity(ind(:,1));
         concat = concat + 1;
        else
         sst = [sst; obs.temperature(ind(:,1))]; sss = [sss; obs.salinity(ind(:,1))];
         lonstats = [lonstats; lonstat(:)]; latstats = [latstats; latstat(:)];
         cruiseDate = [cruiseDate; mean(obs.time)];
        end
        % now start loop to go through each grid point within polygon region
        for mm = 1:length(indgrid)
            usestat = dsearch(lonstat,latstat,TRI,longrid(mm),latgrid(mm));
            TEMP = obs.temperature(ind(usestat,1):ind(usestat,2));
            SALT = obs.salinity(ind(usestat,1):ind(usestat,2));
            DPTH = obs.pressure(ind(usestat,1):ind(usestat,2));
            bad = find(diff(DPTH)==0)+1;
            if(~isempty(bad)); % get rid of duplicates
                TEMP(bad) = []; SALT(bad)=[]; DPTH(bad)=[]; 
            end
            DPTH = sw_dpth(DPTH,mean(latstat));
            % get rid of NaN's
            goodS = find(~isnan(SALT)); goodT = find(~isnan(TEMP));
            SALT=SALT(goodS); DPTH_S = DPTH(goodS);
            TEMP=TEMP(goodT); DPTH_T = DPTH(goodT);
             % pad to 0 depth with shallowest data point
            DPTH_T = [0;DPTH_T];  DPTH_S = [0;DPTH_S];
            SALT = [SALT(1);SALT]; TEMP = [TEMP(1);TEMP];
            % pad to bottom with deepest point
            maxh = h(indgrid(mm)); %depth from grid file
            z_roms = zz(:,indgrid(mm)); % depths of roms z-grid
            AT = interp1(DPTH_T, TEMP, z_roms, 'nearest', 'extrap');
            AS = interp1(DPTH_S, SALT, z_roms, 'nearest', 'extrap');
            [indi,indj] = ind2sub([M L],indgrid(mm));
            salt_roms(:,indi,indj) = AS;
            temp_roms(:,indi,indj) = AT;
            if mm == 500*round(mm/500)
                disp(['      mm = ',num2str(mm)])
            end;
        end
    end
end
cruiseDate(isnan(cruiseDate))=[];
cruiseDate = nanmean(cruiseDate);

if(isempty(cruiseDate))
    disp(['no data available in this time span for any region']);
    cruiseDate = 0;
end

% now write to the ocean_ini file, first getting rid of all NaN's
bad = find(isnan(salt_roms));
% if(~isempty(bad));for i=1:length(bad)
%     salt_roms(bad(i)) = salt_roms(bad(i)-1);
%     end;end
varval(1,:,:,:) = salt_roms; 
nc_varput(outfile, 'salt', varval); clear varval bad
% bad = find(isnan(temp_roms));
% if(~isempty(bad));for i=1:length(bad)
%     temp_roms(bad(i)) = salt_roms(bad(i)-1);
%     end;end
varval(1,:,:,:) = temp_roms; 
nc_varput(outfile, 'temp', varval); clear varval

disp(['...DONE']);

toc

if(cruiseDate~=0) %debugging
%%% then plot to see effects
figure;
  set(gcf,'position',[234 346 1148 603]);
   subplot(121) % plot surface salinity of old IC and stations
    ss = squeeze(salt_old(end,:,:)); ss(~logical(mask_rho))=NaN;
    h1 = pcolorcen(lon_rho,lat_rho,ss); caxis([26 31]);
    hold on; colorbar
    scatter(lonstats,latstats,55,sss,'filled');
    plot_WAcoast;
    axis image; axis([-123.4 -122 47 49.5]);
   subplot(122) % plot surface salinity of new IC and stations
    ss = squeeze(salt_roms(end,:,:)); ss(~logical(mask_rho))=NaN;
    h1 = pcolorcen(lon_rho,lat_rho,ss); caxis([26 31]);
    hold on; colorbar
    scatter(lonstats,latstats,55,sss,'filled');
    plot_WAcoast;
    axis image; axis([-123.4 -122   47  49.5]);
end

    
    
    
    







  
