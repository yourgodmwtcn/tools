function [] = dc_vor_fric(Tdir,infile,basename,tt)
% Plots depth averaged vorticity terms
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready
%%
% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdYlGn',32));
fontSize = [12 10 12];

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

toUTM =  findstr(ncreadatt(infile,'lon_u','units'),'degree');

[vor,xvor,yvor,zvor] = vorticity_cgrid(grid1,ubar,vbar,toUTM);

% total water depth
H = zeta + G.h;
Cd = ncread(infile,'rdrg2');
% centered difference
Hx = diff(H,2,1)./mean(G.DX(:))/2;
Hy = diff(H,2,2)./mean(G.DY(:))/2;

urho = avg1(ubar(:,2:end-1),1);
vrho = avg1(vbar(2:end-1,:),2);
Urho = hypot(urho,vrho);
Hrho = H(2:end-1,2:end-1)';

vor_rho = avg1(avg1(vor,1),2);

% vorticity dissipation
diss  = -Cd*Urho.*vor_rho./Hrho;

% slope torque
slope = -Cd*Urho./(Hrho.^2) ... % Cd |u|/H^2
           .* (urho.*Hy(2:end-1,:)' - vrho.*Hx(:,2:end-1)');% u*Hy - v*Hx
       
% advective term - on PSI points
adv = avg1(ubar,2).*vor + avg1(vbar,1).*vor;

% stretching term
stretch = (vor_rho + 1e-4)./Hrho .* (avg1(ubar(:,2:end-1),1).*Hx(:,2:end-1)' + avg1(vbar(2:end-1,:),2).*Hy(2:end-1,:)');

% speed torque
Upsi = hypot(avg1(ubar,2),avg1(vbar,1)); 
Upsi_x = diff(Upsi,1,1)./mean(G.DX(:));
Upsi_y = diff(Upsi,1,2)./mean(G.DY(:));
speed = Cd.*(avg1(ubar(:,2:end-1).*Upsi_y,1) - avg1(vbar(2:end-1,:).*Upsi_x,2))./Hrho;

%% make plots
h(1) = subplot(241);
imagescnan(xvor,yvor,vor'); shading flat
cvec = [-5 5]*1E-3; caxis(cvec);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Vorticity (s^{-1})','fontweight','bold')
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
imagescnan(G.lon_psi,G.lat_psi,adv'); shading flat
cvec = [-1 1]*1E-3; caxis(cvec);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Advection','fontweight','bold')
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
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),stretch'); shading flat
hold on 
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,[40 80],'k-');
clabel(C,hc,'LabelSpacing',916,'FontWeight','normal','FontSize',10);
cvec = [-5 5]*1e-7; caxis(cvec);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Stretching','fontweight','bold')
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
imagescnan(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),(speed+diss+slope)'); shading flat
cvec = [-5 5]*1e-7;caxis(cvec);
hcbar = colorbar('Southoutside'); 
% fix scaling
Z_dar;
% add labels
title('Friction','fontweight','bold')
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

%suplabel('Depth Averaged Vorticity Equation Terms','t');
% add tides & resize
h(5) = subplot(2,4,[5:8]);

set(h(1),'Units','pixels')
linkprop(h,'Units');

for i=1:4
    set(h(i),'Position',get(h(i),'Position')+[0 -225 0 180]);
end
pos1 = get(h(5),'Position');
set(h(5),'Position',pos1 + [0 0 0 -175])

set(h(1),'Units','normalized')
linkprop(h,'Units');

plot_tides(T.time_datenum); 
beautify(fontSize);
