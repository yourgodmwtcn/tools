function []=energy_flux(M)

% energy_flux.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction

td = M.td;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

z = M.z_rho; zw = M.z_w;
eta = M.zeta; zb = M.z_w(1,:);
u = M.u; v = M.v;
t = M.temp; s = M.salt;

[NR,NC] = size(s);

% CALCULATE IGW ENERGY FLUX
%
% first create the pressure due to density (not due to surface height)
g = 9.8;
z0 = mean(z,2); zw0 = mean(zw,2);
dz = diff(zw0);
dzmat = dz * ones(1,NC);
rho = squeeze(Z_make_potdens(s,t));
rho_fud = flipdim(rho,1);
dzmat_fud = flipdim(dzmat,1);
p_fud = cumsum(g*rho_fud.*dzmat_fud);
p = flipdim(p_fud,1);
p = [p; zeros(1,NC)]; % pressure at interfaces, packed from bottom (z,t)
%
% remove the time average at each depth
pp = p - Z_godin(p')';
%
% interpolate to bin centers
pp_rho = pp(1:end-1,:) + diff(pp,1,1)/2;
%
% create depth average
P = ones(NR,1) * sum(dzmat.*pp_rho,1)/sum(dz);
%
% create deviation from depth average (baroclinic pressure)
pbc = pp_rho - P;
%
% create baroclinic velocity
ubar = ones(NR,1) * M.ubar;
vbar = ones(NR,1) * M.vbar;
ubc = u - ubar;
vbc = v - vbar;
%
% create pressure work
if 1
    % band pass filter first
    M.pbc = pbc;
    M.ubc = ubc;
    M.vbc = vbc;
    FP.vars={'pbc';'ubc';'vbc'};
    M.yday=td0;
    M.z=z0;
    M=Z_FilterCTD2(M,FP);
    pwx = M.pbcf .* M.ubcf;
    pwy = M.pbcf .* M.vbcf;
    tag = 'band-passed';
else
    pwx = pbc .* ubc;
    pwy = pbc .* vbc;
    tag = 'not band-passed';
end
% filter out tides
Fzx = Z_godin(pwx')';
Fzy = Z_godin(pwy')';
%
% take vertical integral
Fx = sum(Fzx.*dzmat,1);
Fy = sum(Fzy.*dzmat,1);

figure; set(gcf,'position',[20 20 1400 900]); Z_fig;

plot(td0,Fx,'-r',td0,Fy,'-b')
aa = axis;
axis([td0(1) td0(end) aa(3) aa(4)]);
[xt,yt] = Z_lab('ul');
text(xt,yt,['(red) Fx = ',num2str(round(nanmean(Fx))),' (W m^{-1}) '])
[xt,yt] = Z_lab('ll');
text(xt,yt,['(blue) Fy = ',num2str(round(nanmean(Fy))),' (W m^{-1}) '])
%
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')
grid on

