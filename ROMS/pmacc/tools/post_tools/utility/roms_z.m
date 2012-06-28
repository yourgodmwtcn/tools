function z = roms_z(H,zeta,cs);

% z = roms_z(H,zeta,cs);
% z = roms_z(H,cs);
%
% returns a 3D z variable given 2D variables H (bathymetry) and zeta (surface height),
% and the vector cs (Cs_r or Cs_r in ROMS terms).
%
% IMPORTANT: note that this only applies for ROMS Vtransform = 1 AND hc = 0
%
% neil banas feb 2009

if nargin == 2
	cs = zeta;
	zeta = zeros(size(H));
end

K = length(cs);
[J,I] = size(H);

H3 = repmat(reshape(H,[1 J I]),[K 1 1]);
zeta3 = repmat(reshape(zeta,[1 J I]),[K 1 1]);
cs3 = repmat(reshape(cs,[K 1 1]),[1 J I]);

z = zeta3 + cs3 .* (H3 + zeta3);
