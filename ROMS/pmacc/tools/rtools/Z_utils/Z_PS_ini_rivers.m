function Z_PS_ini_rivers(year, S, out_dir, ininame, riv_dir, riverFile, riverpolyFile)
% Z_PS_ini_rivers.m DAS, 2010
%
% m-file for adjusting the IC in rivers that were carved out after the IC
% was made originally from the climatology file. It linearly interpolates
% along the river axes in the riverpolyFile to make S go from 0 at head to
% S_ocean at mouth. Temp is interpolated from Temp of river at head to
% T_ocean at mouth.
%
% NEEDS: year, S, out_dir (run directory), ininame (IC file name), riv_dir
% (where rivers.nc is), riverFile (file used to carve channels out),
% riverpolyFile (polygons that enclose each river)
%
% EDIT: 1) do_check == 1 to plot 
%
% written by DAS 3/9/2009 using Prism_Jun_2006.nc as example
% updated by sng 12/21/2011 to be more versatile on file structure
%--------------------------------------------------------------------------

do_check = 0; %plots original IC and altered IC

gfile = [out_dir,'grid.nc']; %assumes already made grid in the run folder
inifile = [out_dir,ininame]; %assumes already made ini in the run folder
tt = nc_varget(inifile,'ocean_time'); 
tt = tt/86400; % Ini date in yeardays (0 based)
h = nc_varget(gfile,'h'); % get grid data 
[M,L] = size(h);
mask_rho = nc_varget(gfile,'mask_rho');
lat_rho = nc_varget(gfile,'lat_rho');
lon_rho = nc_varget(gfile,'lon_rho');

rfile_in = [riv_dir 'PS_flow_',num2str(year) '.mat'];
load(rfile_in);
timeriv = dsearchn(Qr_yearday',tt);

load(riverFile);
num_rivers = length(rivers);
load(riverpolyFile); %needs to have polys struct defined
%polynames = fieldnames(polys); 
%npoly = length(polynames);

%%%% now read in grid info from IC, which is also the outfile
outfile = [out_dir,ininame]; 

% read in initial condition
salt_roms = nc_varget(outfile,'salt'); salt_old =salt_roms; %for debugging plot below
temp_roms = nc_varget(outfile,'temp'); % temp_old =temp_roms; %for debugging plot below

%%% now interpolate from seaward station to innermost station for each river region
 for i=1:num_rivers
    region = rivers(i).name;
   if(~strcmp('skagit_south',region)); %skip skagit south if it's there
    disp([' Working on ', region, ' river IC'])
    xr = rivers(i).lon; yr = rivers(i).lat; %first point is closest to basin
    mr = length(rivers(i).lat); 
    Sal = NaN*ones(S.N,mr); Temp = Sal;
    xx = sw_dist(yr, xr, 'km'); xx=[0;cumsum(xx)];
    % find point in grid closest to xsea, ysea and extract S,T there
    ygrd = dsearchn(lat_rho(:,1),yr(1)); xgrd = dsearchn(lon_rho(1,:)',xr(1)); 
    ss = (nc_varget(outfile,'salt')); 
    Sal(:,1) = squeeze(ss(:,ygrd,xgrd)); Sal(:,end) = 0; 
    tt = (nc_varget(outfile,'temp')); 
    Temp(:,1) = squeeze(tt(:,ygrd,xgrd)); Temp(:,end) = T_riv(i,timeriv);

    % make distance along river channel variable and interpolate from ocean to 0
    for k=1:S.N %interpolate at each N level from ocean value to river value
        Sal(k,2:end-1) = interp1(xx([1 end]),Sal(k,[1 end]),xx(2:end-1));
        Temp(k,2:end-1) = interp1(xx([1 end]),Temp(k,[1 end]),xx(2:end-1));
    end

    % now find points on grid that are within river polygon
    %AA = getfield(polys,region); %lat/lon of river polygon
    AA = polys.(region); % this is the prefered way to get a field using a variable to hold the field name.
    rx = AA(1:end-1,1); ry = AA(1:end-1,2); %exclude last point b/c it's same as first
    inGrid = inpolygon(lon_rho, lat_rho, rx, ry);
    indgrid = find(inGrid==1);
    longrid = lon_rho(inGrid);
    latgrid = lat_rho(inGrid);
    % now do delaunay tri for interpolated data, then loop through
    %TRI = delaunay(xr,yr);
    DT = DelaunayTri(xr,yr);
    for mm = 1:length(indgrid)
            %usestat = dsearch(xr,yr,TRI,longrid(mm),latgrid(mm));
            usestat = nearestNeighbor(DT,longrid(mm),latgrid(mm));
            AAT = Temp(:,usestat);
            AAS = Sal(:,usestat);
            [indi,indj] = ind2sub([M L],indgrid(mm));
            salt_roms(:,indi,indj) = AAS;
            temp_roms(:,indi,indj) = AAT;
            if mm == 50*round(mm/50)
                disp(['      mm = ',num2str(mm)])
            end;
     end
   end %end if skagit south
 end %end for i=1:numrivers

 % now write to the ocean_ini file
varval(1,:,:,:) = salt_roms;
nc_varput(outfile, 'salt', varval); clear varval
varval(1,:,:,:) = temp_roms;
nc_varput(outfile, 'temp', varval); clear varval

disp('...DONE');

if(do_check) %debugging
%%% then plot to see effects
figure;
  set(gcf,'position',[234 346 1148 603]);
   subplot(121) % plot surface salinity of old IC and stations
    ss = squeeze(salt_old(end,:,:)); ss(~logical(mask_rho))=NaN;
    h1 = pcolorcen(lon_rho,lat_rho,ss); caxis([20 31]); %#ok<NASGU>
    hold on; colorbar
    scatter(xr,yr,45,Sal(end,:),'filled');
    plot_WAcoast('detailed');
    axis image; axis([-123.4 -122   47   49.5]);
   subplot(122) % plot surface salinity of new IC and stations
    ss = squeeze(salt_roms(end,:,:)); ss(~logical(mask_rho))=NaN;
    h1 = pcolorcen(lon_rho,lat_rho,ss); caxis([20 31]); %#ok<NASGU>
    hold on; colorbar
    scatter(xr,yr,45,Sal(end,:),'filled');
    plot_WAcoast('detailed');
    axis image; axis([-123.4 -122   47   49.5]);
end

    
    
    
    







  
