function [ape] = ape(rho)
    
    [n,xout] = hist(rho(:),50);
    
    N = numel(rho);
    
    nz = size(rho,2);
    
    nnew = n./N*nz;
    
    nnew .* xout;