% ASSUMES REGULAR GRID
function [floats] = tracmass_read(fname,rgrid)

% from the manual
% The trajectories are stored in <outDataDir>/<outDataFile>_run.asc, which has to be
% specified in <project>_run.in
% The trajectories are stored in colums of
%
%              ntrac,niter,x1,y1,z1,tt,t0,subvol,temp,salt,dens
%
% where
% ntrac is the trajectory number
% niter is the TRACMASS code iteration (only important for TRACMASS modellers)
% x1 is the zoonal position of the trajectory particle - INDEX
% y1 is the meridional position of the trajectory particle - INDEX
% z1 is the vertical position of the trajectory particle - INDEX
% tt the time of the trajectory particle (in days)
% t0 the initial time of the trajectory particle
% subvol is the the "volume" transport in m3/s of the trajectory
% temp is the temperature of the trajectory particle
% salt is the salinity/specific humidity of the trajectory particle
% dens is the density of the trajectory particle
   tic 
   [ntrac,~,ix,iy,iz,tt,t0,subvol,temp,salt,~] = ...
                textread(fname,'%d%d%f%f%f%f%f%f%f%f%s');
   toc
   disp('Finished opening file. Now processing for unique records');
   tic         
   floats.time = unique(tt);
   
   xr = rgrid.x_rho(1,:);
   yr = rgrid.y_rho(:,1)';
   
   dx = xr(2)-xr(1); dy = yr(2)-yr(1);
   
   for i = 1:length(unique(ntrac)) % ith drifter
      k=1;
      for j=1:length(ntrac)    
        if ntrac(j) == i
            floats.t0(i) = t0(j);
            fx = floor(ix(j)); fy = floor(iy(j)); fz = floor(iz(j));
            cz = ceil(iz(j));
            if fy == 0, fy = 1; end
            if fx == 0, fx = 1; end
            dz = rgrid.z_r(cz,fy,fx) - rgrid.z_r(fz,fy,fx);
            % not all floats start at t=0
            dt = floats.t0(i)-floats.time(1);
            floats.x(k+dt,i) = xr(fx) + (ix(j)-fx) * dx;
            floats.y(k+dt,i) = yr(fy) + (iy(j)-fy) * dy;
            floats.z(k+dt,i) = rgrid.z_r(fz,fy,fx) + (iz(j)-fz) * dz;
            floats.temp(k+dt,i) = temp(j);
            floats.salt(k+dt,i) = salt(j);
            floats.t(k+dt,i)    = tt(j);
            floats.subvol(k+dt,i) = subvol(j);
            k=k+1;
        end
      end
   end
   floats.init(:,1) = floats.x(1,:)';
   floats.init(:,2) = floats.y(1,:)';
   floats.init(:,3) = floats.z(1,:)';
   
   % fill with nans
   names = fieldnames(floats);
   for ii=1:length(names)
      floats.(names{ii}) = fillnan(floats.(names{ii}),0); 
   end
   
   floats.fac = 1; % outputs at model output
   toc