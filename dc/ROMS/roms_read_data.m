function [out] = roms_read_data(folder,varname,start,count,stride)
    
    % get all history files
    if isdir(folder)
        files = roms_find_file(folder,'his');
    else
        files = folder;
    end
        
    k = 1;
    
    for ii=1:length(files)
        if isdir(folder)
            fname = [folder '/' char(files(ii))];
        else
            fname = folder;
        end
        if ii == 1
            vinfo  = ncinfo(fname,varname);
            dim = length(vinfo.Size);
            
            if ~exist('start','var'), start = ones([1 dim]); end
            if ~exist('count','var'), count = inf([1 dim]); end
            if ~exist('stride','var'), stride = ones([1 dim]); end
        end
        temp = squeeze(double(ncread(fname,varname,start,count,stride)));
        if count(end) == 1
            out = temp;
            return;
        end
        
        switch ndims(temp)
            case 2
                out(k:k+length(temp)-1) = temp;
                k = k+length(temp);
            case 3
                out(:,:,k:k+size(temp,3)-1) = temp;
                k = k+size(temp,3);
            case 4
                out(:,:,:,k:k+size(temp,4)-1) = temp;
                k = k+size(temp,4);
        end        
    end   
    