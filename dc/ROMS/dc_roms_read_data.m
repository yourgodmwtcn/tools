function [out,xax,yax,zax] = dc_roms_read_data(folder,varname,tindices,volume,stride)
    
    disp(['Reading ' varname]);
    
    % set flags
    k = 1;
    quitflag = 0;
    objflag = 0;
    
    % set inputs
    if ~exist('tindices','var') || isempty(tindices), tindices = [1 Inf]; end
    if ~exist('volume','var') || isempty(volume), volume = {}; end
    if ~exist('stride','var'), stride = [1 1 1 1]; end
    
    if isobject(folder)
        run = folder;
        folder = run.dir;
        objflag = 1;
    end
    
    tic;
    if strcmpi(varname,'vor')
        files{1}='ocean_vor.nc';
    else
        % get all history files
        if isdir(folder)
            files = roms_find_file(folder,'his');
        else
            files = folder;
        end
    end
    
    %unused params
    slab = Inf;
    
    
    for ii=1:length(files)
        if isdir(folder)
            fname = [folder '/' char(files(ii))];
        else
            fname = folder;
        end
        % variable information
        vinfo  = ncinfo(fname,varname);
        nt     = vinfo.Size(4);
        if ii == 1
            dim = length(vinfo.Size);  
            if objflag
                grd = run.rgrid;
            else
                grd = roms_get_grid(fname,fname,0,1);
            end
            [xax,yax,zax,vol] = dc_roms_extract(grd,varname,volume,1);
            %[~,~,~,time,xunits,yunits] = dc_roms_var_grid(grd,varname);
        end
        % process tindices (input) according to file
        [~,tnew,dt,~,tstride] = roms_tindices(tindices,slab,nt);
        stride(4) = tstride(end);
        % Case 1 : if requested data not in this file, skip to next
        if tnew(1) > vinfo.Size(end)
            tindices = tindices - nt;
            continue;
        end
        % Case 2 : requested data spans 2 files
        if tnew(1) < nt && tnew(2) > nt
            % change ending for current read
            tnew(2) = Inf;
            % set next read to start from beginning of new file
            tindices(1) = 1;
        end
        % Case 3 : requested data finishes in current file
        if tnew(2) <= nt && ~isinf(tindices(end))
            quitflag = 1;
        end
        [start,count] = roms_ncread_params(dim,0,1,slab,tnew,dt,vol);
        
        temp = squeeze(double(ncread(fname,varname,start,count,stride)));
        if count(end) == 1
            out = temp;
            return;
        end
        
        if isvector(temp)
            out(k:k+length(temp)-1) = temp;
            k = k+length(temp);
        else        
            switch ndims(temp)
                case 2
                    out(:,k:k+size(temp,2)-1) = temp;
                    k = k+size(temp,2);
                case 3
                    out(:,:,k:k+size(temp,3)-1) = temp;
                    k = k+size(temp,3);
                case 4
                    out(:,:,:,k:k+size(temp,4)-1) = temp;
                    k = k+size(temp,4);
            end  
        end
        if quitflag, break; end
    end 
    toc;
    