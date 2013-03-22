function [out] = roms_read_data(folder,varname,tindices)
    
    % get all history files
    folder = 'runs\runte-11\';
    files = ls([folder '*_his*.nc']);
        
    k = 1;
    
    for ii=1:size(files,1)
        fname = [folder '\' files(ii,:)];
        
        temp = double(ncread(fname,varname));
        if ndims(temp) == 3
            out(:,:,k:k+size(temp,3)-1) = temp;
            k = k+size(temp,3);
        else
            out(:,:,k:k+size(temp,4)-1) = temp;
            k = k+size(temp,4);
        end
        
    end   
    