function [] = dc_bernoulli_velvec(Tdir,infile,basename,tt)
% <Modified from dc_bernoulli.m
% Modified Parker's basic.m code
% 10 July 2012 - Deepak Cherian

% get file info
[G,S,T] = Z_get_basic_info(infile);

cmap = flipud(lbmap(32,'RedBlue'));
cmap2 = cmap;
fontSize = [18 18 24];

lon0 = 0;-122.95;
lat0 = 0;48.465; 

% from http://www.csgnetwork.com/degreelenllavcalc.html
lat2km = 1;111199.32/1000;
lon2km = 1;73952.24/1000;

%% plot gz
clf
cvec = [0 8]; 
zt = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);
zt = zt - nanmean(zt(:));
zt = 9.81*(zt - nanmin(zt(:)));

% plot 1/2 (u^2+v^2) + g zeta
u = double((ncread(infile,'u',[1 1 S.N 1],[Inf Inf 1 Inf])));
v = double((ncread(infile,'v',[1 1 S.N 1],[Inf Inf 1 Inf])));
uavg = avgx(u.*u);
vavg = avgy(v.*v);
u2v2 = squeeze(0.5*(uavg(:,2:end-1) + vavg(2:end-1,:)));

%h(2) = subplot(232);
pcolorcen(G.lon_rho(2:end-1,2:end-1)-lon0,G.lat_rho(2:end-1,2:end-1)-lat0,(u2v2+zt(2:end-1,2:end-1)')'); 
shading flat
colormap(cmap2);
caxis(cvec);
hcbar = colorbar('EastOutside'); 
% fix scaling
Z_dar;
% add time info
%Z_info(basename,tt,T,'lr'); 
% add labels
titlestr = sprintf('\\emph{B} (m$^2$/s$^2$) : t =  %.2f flood', ...
    (T.time_datenum -  7.328527395833334e+05)/datenum(0,0,0,7,0,0));
ht = title(titlestr,'fontweight','bold');
set(ht,'interpreter','latex');
% add coastline
Z_addcoast('detailed',Tdir.coast);
hold on
% add transect lines
transectx = ([-122-56/60-49/3600 -122-59/60-20/3600; ...
             -122-56/60-49/3600 -122-59/60-49/3600; ...
             -122-56/60-41/3600 -123-00/60-07/3600] - lon0)*lon2km;

transecty = ([48+28/60+22/3600 48+28/60+22/3600; ...
             48+28/60+51/3600 48+28/60+51/3600; ...
             48+29/60+25/3600 48+29/60+25/3600] - lat0)*lon2km;
fmt = 'k-';
plot(transectx(1,:),transecty(1,:),fmt,'LineWidth',3);
plot(transectx(2,:),transecty(2,:),fmt,'LineWidth',3);
plot(transectx(3,:),transecty(3,:),fmt,'LineWidth',3);

% and velocity vectors
Z_velvec(infile,G,S,'mr')
xlabel('Longitude (deg)'); ylabel('Latitude (deg)');
ax(2) = gca; 
set(gca,'XTick',sort([-123.03:0.03:128.96]));
beautify(fontSize); box on
set(gcf,'renderer','zbuffer');

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;
