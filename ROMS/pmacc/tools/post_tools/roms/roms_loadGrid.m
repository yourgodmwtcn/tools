function G = roms_loadGrid(A);

% G = roms_loadGrid(filename);
% G = roms_loadGrid(seriesDef);
%
% neil banas mar 2009

if ischar(A)
	fname = A;
elseif isfield(A,'dirname') & isfield(A,'basename') & isfield(A,'ncn')
	fname = roms_filename([A.dirname A.basename],A.ncn(1));
else
	error('roms_loadGrid requires either a filename or a seriesDef.');
end

G.lon = nc_varget(fname, 'lon_rho');
G.lat = nc_varget(fname, 'lat_rho');
% PM additions 3/7/2012
pm = nc_varget(fname, 'pm');
pn = nc_varget(fname, 'pn');
G.DX = 1./pm;
G.DY = 1./pn;
% end PM edit
G.lonu = nc_varget(fname, 'lon_u');
G.latu = nc_varget(fname, 'lat_u');
G.lonv = nc_varget(fname, 'lon_v');
G.latv = nc_varget(fname, 'lat_v');
try
	G.cs = nc_varget(fname, 'Cs_r');
	G.csw = nc_varget(fname, 'Cs_w');
end

G.mask = nc_varget(fname,'mask_rho');
G.masku = nc_varget(fname,'mask_u');
G.maskv = nc_varget(fname,'mask_v');

G.H = nc_varget(fname, 'h');
G.Hu = interp2(G.lon,G.lat,G.H,G.lonu,G.latu);
G.Hv = interp2(G.lon,G.lat,G.H,G.lonv,G.latv);
%G.H(G.mask==0) = nan; %change by DAS because getting NaN's in points that were too close to grid boundary
%G.Hu(G.masku==0) = nan;
%G.Hv(G.maskv==0) = nan;


% K,J,I in Neil's jargon = N,L,M in ROMS jargon
%G.K = length(G.cs);
%[G.J, G.I] = size(G.lon);