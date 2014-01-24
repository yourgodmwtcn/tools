%function [] = fancy(Tdir,infile,basename,tt);

clear
close all
Tdir = pan_start;
infile = ['/Users/PM2/Documents/salish_2006_4/OUT', ...
    '/ocean_his_5071.nc'];

%
% plots Salish simulations, a fancy plot
% 8/8/2012  Parker MacCready

% get file info
[G,S,T]=Z_get_basic_info(infile);

% plot salinity
s0 = nc_varget(infile,'salt',[0 S.N-1 0 0],[1 1 -1 -1]);
% eta = nc_varget(infile,'zeta');
% surface(G.lon_rho,G.lat_rho,eta,s0);
% shading interp
% cvec = [26 32]; caxis(cvec);

surfl(G.lon_rho,G.lat_rho,s0)
shading interp
colormap copper

aa = axis;
meanlat = mean(aa(3:4));
dar = [1/cos(pi*meanlat/180) 1 100]; % Cartesian scaling
set(gca,'dataaspectratio',dar);

view(3);