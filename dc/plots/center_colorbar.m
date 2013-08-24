function [] = center_colorbar()

    colorbar;
    clim = caxis;
    
    a = max(abs(clim));
    caxis([-a a]);
    