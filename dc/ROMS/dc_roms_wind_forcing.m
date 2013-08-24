function [] = dc_roms_wind_forcing(S,frcname)

    if exist(frcname,'file'), delete(frcname); end

    nccreate(frcname,'sms_time','Dimensions',{'sms_time' Inf});
    nccreate(frcname,'sustr','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2 'sms_time'});
    nccreate(frcname,'svstr','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2 'sms_time'});
    
    nccreate(frcname,'xi_rho','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2'});
    nccreate(frcname,'eta_rho','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2});
    
    ncwrite(frcname,'xi_rho',S.x_rho);
    ncwrite(frcname,'eta_rho',S.y_rho);
    
    % add attributes
    ncwriteatt(frcname,'sustr','long_name','surface u-momentum stress');
    ncwriteatt(frcname,'sustr','units', 'Newton meter-2');
    ncwriteatt(frcname,'sustr','coordinates', 'xi_rho eta_rho sms_time');
    ncwriteatt(frcname,'sustr','time','sms_time');
    
    ncwriteatt(frcname,'svstr','long_name','surface v-momentum stress');
    ncwriteatt(frcname,'svstr','units', 'Newton meter-2');
    ncwriteatt(frcname,'svstr','coordinates', 'xi_rho eta_rho sms_time');
    ncwriteatt(frcname,'svstr','time','sms_time');
    
    ncwriteatt(frcname,'sms_time','long_name','seconds since 0001-01-01 00:00:00');
    ncwriteatt(frcname,'sms_time','units','seconds');
    
    ncwriteatt(frcname,'xi_rho','long_name','X-location of RHO-points');
    ncwriteatt(frcname,'xi_rho','units','meter');
    ncwriteatt(frcname,'eta_rho','long_name','Y-location of RHO-points');
    ncwriteatt(frcname,'eta_rho','units','meter');
    