function [floats] = roms_read_floats(flt)

    floats.x = ncread(flt,'x')';
    floats.y = ncread(flt,'y')';
    floats.z = ncread(flt,'depth')';
    floats.time = ncread(flt,'ocean_time');
    floats.temp = ncread(flt,'temp')';
    floats.salt = ncread(flt,'salt')';
    
    for ii = 1:size(floats.x,2)
        clear ind
        ind = find(floats.x(:,ii) > 0);
        ind = ind(1);
        
        floats.init(ii,:) = [floats.x(ind,ii) floats.y(ind,ii) floats.z(ind,ii) floats.time(ind)];
    end
    
    floats.comment = 'init = (x,y,z,t) = initial location, release time in meters, seconds';