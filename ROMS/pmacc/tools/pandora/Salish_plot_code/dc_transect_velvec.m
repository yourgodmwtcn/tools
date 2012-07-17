function [] = dc_mom_full(Tdir,infile,basename,tt)
    
% Plots transect velocity vectors and full field - Deepak Cherian
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
%t0 = nc_varget(infile,'temp',[0 S.N-1 0 0],[1 1 -1 -1]);
ubar      = double(squeeze(ncread(infile,'ubar',[1 1 1],[Inf Inf Inf])));
vbar      = double(squeeze(ncread(infile,'vbar',[1 1 1],[Inf Inf Inf])));

grid1.xv = grid.lon_v(1,:)';
grid1.yv = grid.lat_v(:,1);
grid1.zv = grid.z_v(end,1,1);

grid1.xu = grid.lon_u(1,:)';
grid1.yu = grid.lat_u(:,1);
grid1.zu = grid.z_u(end,1,1);

toUTM =  findstr(ncreadatt(infile,'lon_u','units'),'degree');

[vor,xvor,yvor,zvor] = vorticity_cgrid(grid1,ubar,vbar,toUTM);

transects = [-122-56/60-49/3600 -122-59/60-20/3600 48+28/60+22/3600 ; ...
             -122-56/60-49/3600 -122-59/60-49/3600 48+28/60+51/3600 ; ...
             -122-56/60-41/3600 -123-00/60-07/3600 48+29/60+25/3600 ];
         
transectx = [-122-56/60-49/3600 -122-59/60-20/3600; ...
             -122-56/60-49/3600 -122-59/60-49/3600; ...
             -122-56/60-41/3600 -123-00/60-07/3600];

transecty = [48+28/60+22/3600 48+28/60+22/3600; ...
             48+28/60+51/3600 48+28/60+51/3600; ...
             48+29/60+25/3600 48+29/60+25/3600];
         
x1 = transects(1,1):-0.001:transects(1,2); y1 = repmat(transects(1,3),size(x1));
x2 = transects(2,1):-0.001:transects(2,2); y2 = repmat(transects(2,3),size(x2));
x3 = transects(3,1):-0.001:transects(3,2); y3 = repmat(transects(3,3),size(x3));

%% make plots
clf
% plot full velocity field
h(1) = subplot(221);
imagescnan(xvor,yvor,vor'); shading flat
cvec = [-5 5]*1E-3; caxis(cvec);
hcbar = colorbar('Eastoutside'); 
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
plot(transectx(1,:),transecty(1,:),'k-','LineWidth',1.5);
plot(transectx(2,:),transecty(2,:),'k-','LineWidth',1.5);
plot(transectx(3,:),transecty(3,:),'k-','LineWidth',1.5);
beautify(fontSize); box on;

aa = axis;darscale = 1/cos(pi*mean(aa(3:4))/180); Dlat = aa(4)-aa(3); Dlon = aa(2)-aa(1);
ufact = Dlat/30;
% plot transects
h(2) = subplot(222);
[C,hc] = contour(G.lon_rho,G.lat_rho,G.h,'k');hold on;
ylim([48.44 48.5])
clabel(C,hc,'LabelSpacing',720,'FontWeight','normal','FontSize',10);
quiver(x1,y1,interp2(G.lon_u,G.lat_u,ubar'*darscale*ufact,x1,y1),interp2(G.lon_v,G.lat_v,vbar'*ufact,x1,y1),0,'.','LineWidth',1.5,'Color',[0    0.4549    0.7373]);
quiver(x2,y2,interp2(G.lon_u,G.lat_u,ubar'*darscale*ufact,x2,y2),interp2(G.lon_v,G.lat_v,vbar'*ufact,x2,y2),0,'.','LineWidth',1.5,'Color',[0    0.4549    0.7373]);
quiver(x3,y3,interp2(G.lon_u,G.lat_u,ubar'*darscale*ufact,x3,y3),interp2(G.lon_v,G.lat_v,vbar'*ufact,x3,y3),0,'.','LineWidth',1.5,'Color',[0    0.4549    0.7373]);
[xt,yt] = Z_lab('mr'); uscale = 1;
        deltax = -Dlon/4; deltay = Dlat/20;
        quiver(xt+deltax,yt+deltay,ufact*uscale*darscale,ufact*0,0,'k');
        quiver(xt+deltax,yt+deltay,ufact*0*darscale,ufact*uscale,0,'k');
        text(xt+deltax,yt,[num2str(uscale),' m s^{-1}'],'color','k');
Z_dar;
xlabel('Longitude (deg)')
ylabel('Latitude (deg)');
title('Transects in model');
% add coastline
[handle] = Z_addcoast('detailed',Tdir.coast);
set(handle,'LineWidth',1.5);
beautify(fontSize); box on;

% add tides & resize
h(3) = subplot(2,2,[3 4]);
set(h(1),'Units','pixels')
set(h(1),'Position',get(h(1),'Position')+[0 -225 0 180]);
set(h(2),'Units','pixels')
set(h(2),'Position',get(h(2),'Position')+[0 -225 0 180]);
set(h(3),'Units','pixels')
pos1 = get(h(3),'Position');
set(h(3),'Position',pos1 + [0 0 0 -175])

set(h(1),'Units','normalized')
set(h(2),'Units','normalized')
set(h(3),'Units','normalized')

plot_tides(T.time_datenum); 
beautify(fontSize);

set(gcf,'renderer','zbuffer');


