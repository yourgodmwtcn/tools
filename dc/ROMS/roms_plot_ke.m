function [KE] = roms_plot_ke(fname)

    u = double(ncread(fname,'u'));
    v = double(ncread(fname,'v'));
    
    
    u = avg1(u(:,2:end-1,:,:),1);
    v = avg1(v(2:end-1,:,:,:),2);
    
    [xax,yax,zax,tax,~,~] = dc_roms_var_grid(fname,'rho');
    xax = xax(2:end-1,1,1);
    yax = yax(1,2:end-1,1)';
    
    zax = zax(2:end-1,2:end-1,:);
    
    KE = domain_integrate(0.5 *(u.^2 + v.^2),xax,yax,zax);
    
    plot(tax/86400,KE);
    xlabel('time (days)');
    ylabel('KE / unit vol.');