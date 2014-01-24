function [] = dc_stream_mom(Tdir,infile,basename,tt)
    
% Plots stream-normal depth averaged momentum balance - Deepak Cherian3
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready
%%
% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = flipud(lbmap(32,'RedBlue'));
fontSize = [12 10 14];

% plot vorticity terms
warning off
grid = roms_get_grid(infile,infile,0,1);
warning on

ubar = double(squeeze(ncread(infile,'ubar',[1 1 1],[Inf Inf Inf])));
vbar = double(squeeze(ncread(infile,'vbar',[1 1 1],[Inf Inf Inf])));
zeta = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);

% total water depth
H = zeta + G.h;
Cd = ncread(infile,'rdrg2');
g = 9.81;

urho = avg1(ubar(:,2:end-1),1);
vrho = avg1(vbar(2:end-1,:),2);
Urho = hypot(urho,vrho);
Hrho = H(2:end-1,2:end-1)';

% advective term - on PSI points
uadv = urho .* diff(ubar(:,2:end-1),1,1)./mean(G.DX(:)) ...
        + vrho .* avg1(avg1(diff(ubar,1,2),1),2)./mean(G.DY(:)); % u*u_x + v*u_y

vadv = urho .* avg1(avg1(diff(vbar,1,1),1),2)./mean(G.DX(:)) ...
        + vrho .* diff(vbar(2:end-1,:),1,2)./mean(G.DY(:)); 

% Bottom Friction
u_bfric = -Cd*urho.*Urho./Hrho;
v_bfric = -Cd*vrho.*Urho./Hrho;

% SSH slope
gradPx = -g*avg1(diff(zeta',1,1),1)./mean(G.DX(:));
gradPy = -g*avg1(diff(zeta',1,2),2)./mean(G.DY(:));

alpha = angle(urho + 1i*vrho);
% transform 
[us,~]      = transform(urho,vrho,alpha);
[sadv,nadv] = transform(uadv,vadv,alpha);
[sbfric,~]  = transform(u_bfric,v_bfric,alpha);
[dpds,dpdn] = transform(gradPx(:,2:end-1),gradPy(2:end-1,:),alpha);

%% make plots

clf

cvec_u = [-2 2]; % for velocity field
cvec_ssh = [-0.35 0.35];
cvec_acc = [-2 2]*1e-3; % for acceleration terms
cvec_acc2 = [-5 5]*1e-4;

h(1) = subplot(241);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),us'); shading flat
caxis(cvec_u);
%hcbar = colorbar('Westoutside'); 
% fix scaling
Z_dar;
% add labels
title('u_s')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)');
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'mr')
colormap(cmap);
lx = xlim;
ax(1) = gca; set(gca,'XTick',[lx(1):0.05:lx(2)]);
beautify(fontSize); box on;

h(2) = subplot(242);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),sadv'); shading flat
caxis(cvec_acc);
%hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Along-stream Advection','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
set(gca,'YTickLabel',[]);
ax(2) = gca;
colormap(cmap);
beautify(fontSize); box on;

h(3) = subplot(243);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),dpds'); shading flat
caxis(cvec_acc);
% fix scaling
Z_dar;
% add labels
title('dp/ds','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
set(gca,'YTickLabel',[]);
colormap(cmap);
ax(3) = gca;
beautify(fontSize); box on;

h(4) = subplot(244);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),sbfric'); shading flat
caxis(cvec_acc2);
%hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
hold on
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,[40 80],'k-');
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
% add labels
title('Along stream Friction','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
set(gca,'YAxisLocation','Right');
colormap(cmap);
ax(4) = gca;
beautify(fontSize); box on;

linkaxes(ax,'xy');
linkprop(ax,'xtick');

% make plots - others
h(5) = subplot(245);
imagescnan(G.lon_rho,G.lat_rho,zeta-nanmean(zeta(:))); shading flat
caxis(cvec_ssh);
%hcbar = colorbar('Westoutside'); 
% fix scaling
Z_dar;
% add labels
title('SSHA (m)','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)');
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
colormap(cmap);
lx = xlim;
ax(1) = gca; set(gca,'XTick',[lx(1):0.05:lx(2)]);
beautify(fontSize); box on;

h(6) = subplot(246);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),nadv'); shading flat
caxis(cvec_acc);
%hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar; hold on
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,[40 80],'k-');
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
% add labels
title('Centrifugal','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
set(gca,'YTickLabel',[]);
ax(2) = gca;
colormap(cmap);
beautify(fontSize); box on;

h(7) = subplot(247);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),dpdn'); shading flat
caxis(cvec_acc);
% fix scaling
Z_dar;
% add labels
title('dp/dn','FontWeight','Bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
set(gca,'YTickLabel',[]);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
colormap(cmap);
ax(4) = gca;
beautify(fontSize); box on;

h(8) = subplot(248);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),1e-4*us'); shading flat
caxis(cvec_acc2);
hold on 
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
%hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Coriolis','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
set(gca,'YAxisLocation','Right');
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
colormap(cmap);
ax(3) = gca;
Z_info(basename,tt,T,'lr'); 
beautify(fontSize); box on;

% u_s
handle1 = subplot_cbar(cvec_u,[0.001 0.58 0.09 0.35]);
% zeta
handle2 = subplot_cbar(cvec_ssh,[0.001 0.11 0.09 0.35]);
% adv/dp terms
handle3 = subplot_cbar(cvec_acc,[0.23 0.05 0.1 0.9]);
% friction/coriolis terms
handle3 = subplot_cbar(cvec_acc2,[0.64 0.05 0.1 0.9]);

set(handle1,'YAxisLocation','left')
linkprop([handle1 handle2],'YAxisLocation');

linkaxes(ax,'xy');
linkprop(ax,'xtick');

function [trans_s,trans_n] = transform(inx,iny,alpha)
    % transform to along stream and cross stream co-ordinates
    % eqns. A3 and A4 in Hench and Luettich (2003)
    trans_s = inx.*cos(alpha) + iny.*sin(alpha);
    trans_n = iny.*cos(alpha) - inx.*sin(alpha);
    
