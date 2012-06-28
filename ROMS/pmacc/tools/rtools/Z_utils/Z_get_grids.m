function gridinfo = Z_get_grids(g_dir)
%-----------------------------------------------------
% gridinfo = Z_get_grids(g_dir)
%
% replaces lines of code in make_a_run.m that get available grids in g_dir
% and their properties, listing out on the screen
%
% DAS, Oct 2010
%-----------------------------------------------------

A = dir([g_dir,'*.nc']);

disp('****** Available Grids *************')

for ii = 1:length(A)
    A(ii).sph = nc_varget([g_dir,A(ii).name],'spherical');
    disp([' ',num2str(ii),' = ',A(ii).name,', spherical = ',A(ii).sph]);
end

ng = input('Which grid? (e.g. 1) ');

if ng > length(A); disp('ERROR: out of range'); return; end;

g_name = A(ng).name;
gg_name = strrep(g_name,'.nc','');
sph = A(ng).sph;

% make output
gridinfo.sph = sph;
gridinfo.gfile = [g_dir,g_name];
gridinfo.gg_name = gg_name; 
gridinfo.A = A;