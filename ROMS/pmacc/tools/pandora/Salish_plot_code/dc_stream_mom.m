function [] = dc_stream_mom(Tdir,infile,basename,tt)
    
% Plots stream-normal depth averaged momentum balance - Deepak Cherian3
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready
%%
% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdYlGn',32));
fontSize = [12 10 14];

% plot vorticity terms
warning off
grid = roms_get_grid(infile,infile,0,1);
warning on
%t0 = nc_varget(infile,'temp',[0 S.N-1 0 0],[1 1 -1 -1]);
ubar      = double(squeeze(ncread(infile,'ubar',[1 1 1],[Inf Inf Inf])));
vbar      = double(squeeze(ncread(infile,'vbar',[1 1 1],[Inf Inf Inf])));
zeta      = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);

grid1.xv = grid.lon_v(1,:)';
grid1.yv = grid.lat_v(:,1);
grid1.zv = grid.z_v(end,1,1);

grid1.xu = grid.lon_u(1,:)';
grid1.yu = grid.lat_u(:,1);
grid1.zu = grid.z_u(end,1,1);

% total water depth
H = zeta + G.h;
Cd = ncread(infile,'rdrg2');
g = 9.81;
% centered difference
%Hx = diff(H,2,1)./mean(G.DX(:))/2;
%Hy = diff(H,2,2)./mean(G.DY(:))/2;

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
gradPx = -g*diff(zeta,2,1)./mean(G.DX(:))/2;
gradPy = -g*diff(zeta,2,2)./mean(G.DX(:))/2;

alpha = angle(urho + 1i*vrho);
% transform 
[us,~]     = transform(urho,vrho,alpha);
[sadv,nadv] = transform(uadv,vadv,alpha);
[sbfric,~] = transform(u_bfric,v_bfric,alpha);
[dpds,dpdn] = transform(gradPx(:,2:end-1)',gradPy(2:end-1,:)',alpha);

%% make plots - s

clf

cvec_u = [-2 2]; % for velocity field
cvec_acc = [-5 5]*1e-4; % for acceleration terms
cvec_acc2 = [-5 5]*1e-4;

h(1) = subplot(241);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),us'); shading flat
caxis(cvec_u);
hcbar = colorbar('Southoutside'); 
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
hcbar = colorbar('Southoutside'); 
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
caxis(cvec_acc2);
hold on 
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,[40 80],'k-');
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('dp/ds','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
set(gca,'YTickLabel',[]);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
colormap(cmap);
ax(3) = gca;
beautify(fontSize); box on;

h(4) = subplot(244);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),sbfric'); shading flat
caxis(cvec_acc2);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
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
Z_info(basename,tt,T,'lr'); 
beautify(fontSize); box on;

linkaxes(ax,'xy');
linkprop(ax,'xtick');

% make plots - others
h(5) = subplot(245);
imagescnan(G.lon_v,G.lat_v,vbar'); shading flat
caxis(cvec_u);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('vbar','fontweight','bold')
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

h(6) = subplot(246);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),nadv'); shading flat
caxis(cvec_acc);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Centrifugal (nadv)','fontweight','bold')
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
caxis(cvec_acc2);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('dp/dn','FontWeight','Bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
set(gca,'YAxisLocation','Right');
colormap(cmap);
ax(4) = gca;
Z_info(basename,tt,T,'lr'); 
beautify(fontSize); box on;

h(8) = subplot(248);
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),1e-4*us'); shading flat
caxis(cvec_acc2);
hold on 
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,[40 80],'k-');
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Coriolis','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
set(gca,'YTickLabel',[]);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
colormap(cmap);
ax(3) = gca;
beautify(fontSize); box on;

linkaxes(ax,'xy');
linkprop(ax,'xtick');

function [trans_s,trans_n] = transform(inx,iny,alpha)
    % transform to along stream and cross stream co-ordinates
    % eqns. A3 and A4 in Hench and Luettich (2003)
    trans_s = inx.*cos(alpha) + iny.*sin(alpha);
    trans_n = iny.*cos(alpha) - inx.*sin(alpha);
    
