function [Sint,Savg] = roms_depthIntegrate(S, cs, csw, H, zeta, depthRange);

% [Sint,Savg] = roms_depthIntegrate(S, cs, csw, H, zeta, depthRange);
%
% vertically integrates the variable S between depthRange = [minDepth maxDepth].
% S can be on either the rho grid or the w grid in the vertical. If zeta isn't
% given, uses 0.
% both cs and csw (see roms_loadGrid.m) are required. H and zeta should match S
% in horizontal dimensions but should not be replicated in the vertical dimension.
%
% returns both the depth integral and the depth average.
%
% neil banas mar 2009

S = squeeze(S);
K = size(S,1);
if K==length(cs)
	cstop = csw(2:end);
	csbot = csw(1:end-1);
elseif K==length(csw)
	cstop = [cs 0];
	csbot = [-1 cs];
else
	error('bad size');
end
cstop3 = repmat(cstop(:),[1 size(H)]);
csbot3 = repmat(csbot(:),[1 size(H)]);
zeta3 = repmat(reshape(zeta,[1 size(zeta)]),[K 1 1]);
H3 = repmat(reshape(H,[1 size(H)]),[K 1 1]);
ztop = zeta3 + cstop3.*(zeta3+H3);
zbot = zeta3 + csbot3.*(zeta3+H3);

zmin = min(depthRange) + zeta3;
zmax = max(depthRange) + zeta3;

dz = zeros(size(zmin));
inrange = (zmin <= ztop & zmax >= zbot);
dz(inrange) = min(zmax(inrange),ztop(inrange)) - max(zmin(inrange), zbot(inrange));

S(~isfinite(S)) = 0;
Sint = squeeze(sum(S.*dz));
Savg = Sint ./ squeeze(sum(dz));
Savg(~isfinite(Savg)) = nan;
