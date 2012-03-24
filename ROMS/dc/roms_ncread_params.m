
function [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,vol)

    read_start = ones(1,dim);
    read_count = Inf(1,dim);

    if i == (iend-1)
        read_count(end) = floor((tindices(2)-slab*(i))/dt);
    else
        read_count(end) = floor(slab/dt);%ceil(slab*(i+1)/dt);
    end

    if i == 0
        read_start(end) = tindices(1);
    else
        read_start(end) = slab*i + 1;
    end

    if (iend-1) == 0, read_count(end) = floor((tindices(2)-tindices(1))/dt)+1; end 
    
    % parse vol
    if exist('vol','var')
        read_start(1:end-1) = vol(1:dim-1,1)';
        read_count(1:end-1) = vol(1:dim-1,2)' - vol(1:dim-1,1)' + 1;
    end
