function [out] = roms_read_data(folder,varname,tindices)
    
    % get all history files
    files = ls([folder '*_his*.nc']);
    if isempty(files), files = ls([folder '*_avg*.nc']); end
        
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
    