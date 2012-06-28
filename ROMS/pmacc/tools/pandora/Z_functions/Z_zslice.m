function [interpmat] = Z_zslice(infile,zlev,whichgrid)
% Z_zslice.m 3/2/2011 Parker MacCready

% this returns an array with information to interpolate a ROMS
% file onto a specified z level

% whichgrid can be 'rho' 'u' or 'v'

% get basic grid info
[G,S,T] = Z_get_basic_info(infile);

switch whichgrid
    case 'rho'
        h = G.h;
        mask = G.mask_rho;
    case 'u'
        h = interp2(G.lon_rho,G.lat_rho,G.h,G.lon_u,G.lat_u);
        mask = G.mask_u;
    case 'v'
        h = interp2(G.lon_rho,G.lat_rho,G.h,G.lon_v,G.lat_v);
        mask = G.mask_v;
end
        
[M,L] = size(h);

z_rho = roms_z(h,0*h,S.Cs_r);

zbot = -h;
zbot(~mask)=10;

interpmat = zeros(S.N+2,M,L);

for jj = 1:M
    for ii = 1:L
        % check to see if we should do a calculation
        this_zbot = zbot(jj,ii);
        this_z = squeeze(z_rho(:,jj,ii));
        if this_zbot<zlev
            % add top and bottom
            this_z = [this_zbot; this_z; 0];
            klo = find(this_z<=zlev); klo = klo(end);
            khi = klo+1;
            zlo = this_z(klo); zhi = this_z(khi);
            dz = zhi-zlo;
            fraclo = (zhi-zlev)/dz;
            frachi = 1-fraclo;
            interpmat(klo,jj,ii) = fraclo;
            interpmat(khi,jj,ii) = frachi;
        else
            interpmat(:,jj,ii) = NaN;
        end
    end
end

