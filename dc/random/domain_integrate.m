function [out] = domain_integrate(in,xax,yax,zax)

    out = squeeze(trapz(xax,trapz(yax,trapz(zax,in,3),2),1));