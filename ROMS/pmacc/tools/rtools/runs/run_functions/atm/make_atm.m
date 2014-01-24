function [] = make_atm(gridfile,tdlims,indir0,outdir)
% 1/2/2013  Parker MacCready

year0str = datestr(tdlims(1),'yyyy')
year0 = str2num(year0str);

% make the output directory if needed
if ~exist(outdir,'dir'); mkdir(outdir); end;

% list of input variables
varlist = {'psfc';'t2';'qair';'rain';'swdown';'lwdown';'u10r';'v10r'};

debug = 1;
if debug; varlist = {'psfc','t2','swdown'}; end;

for vv = 1:length(varlist);
    
    tic
    var = varlist{vv};
    
    disp(['* Writing forcing file for ',var]);
    % *** create and populate the NetCDF file for ROMS ***
    
    % get the time interval for this variable (assumes January exists!)
    indir = [indir0,year0str,'/'];
    month = 1;
    load([indir,var,'_',num2str(month),'.mat'],'TD');
    DTD = TD(2)-TD(1);

    % here is vector of times we want, padded one DTD beyond the
    % year ends, it can be different for each variable
    tdvec = [tdlims(1)-DTD:DTD: ...
        tdlims(2)+DTD]';
    
    if debug; tdvec = [tdvec(1:10); tdvec(end-10:end)]; end;
    
    % create index vectors to help find things
    % (note that the processed atm files are stored in month files)
    %
    td = []; % initialize datenum vector
    mo = []; % initialize a month vector
    yr = []; % initialize a year vector
    ind = []; % index within the month
    %
    for year = year0-1:year0+1
        indir = [indir0,num2str(year),'/'];
        if exist(indir,'dir')==7
            for month = 1:12
                infile = [indir,var,'_',num2str(month),'.mat'];
                if exist(infile)==2
                    load(infile,'TD');
                    nt = length(TD);
                    td = [td; TD];
                    mo = [mo; month*ones(nt,1)];
                    yr = [yr; year*ones(nt,1)];
                    ind = [ind; [1:nt]'];
                end
            end % end of month loop
        end
    end % end of year loop
    
    % Read in the grid
    lon_rho = nc_varget(gridfile, 'lon_rho');
    lat_rho = nc_varget(gridfile, 'lat_rho');
    [M,L] = size(lon_rho);
    
    % vartime units = seconds from the start of the year
    vartime = round(86400*(tdvec - datenum(year0,1,1,0,0,0)));
    
    [ncvarname,nclongname,ncunits,nctimename, ...
        scalefactor,scalefactor2] = atm_attributes(var);
    
    % create the NetCDF file
    frcname = [outdir, ncvarname, '.nc'];
    
    disp('  ++ Creating the NetCDF file')
    my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
    nc_create_empty( frcname, my_mode );
    % global attributes
    nc_padheader (frcname, 20000 );
    nc_attput(frcname, nc_global, 'type','ROMS Surface Forcing File');
    
    % define dimensions
    nc_add_dimension(frcname, 'xi_rho', L);
    nc_add_dimension(frcname, 'eta_rho', M);
    nc_add_dimension(frcname, nctimename, length(vartime));
    
    % define variables and attributes
    %   first time
    varstruct.Name = nctimename;
    varstruct.Dimension = {nctimename};
    long_name = [nclongname,' time'];
    units = 'seconds';
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    %   then the field
    varstruct.Name = ncvarname;
    varstruct.Dimension = {nctimename, 'eta_rho', 'xi_rho'};
    long_name = nclongname;
    units = ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    % write the variables
    disp('  ++ Writing to the file')
    nc_varput(frcname, nctimename, vartime);
    
    %BEGIN TIME LOOP
    infile_prev = ['blank'];
    for iii = 1:length(vartime)
        
        % specify the input file
        this_td = tdvec(iii);
        
        % this call to dsearchn gets the closest time
        itd = dsearchn(td,this_td);
        
        if debug
            disp(' ')
            disp(['TARGET: this_td = ',datestr(this_td,0)]);
            disp(['ACTUAL: td(itd) = ',datestr(td(itd),0)]);
        end

        indir = [indir0,num2str(yr(itd)),'/'];
        
        infile = [indir,var,'_',num2str(mo(itd)),'.mat'];
        
        % only reload if the name of indir changes
        if ~strcmp(infile,infile_prev)
            disp(['   loading (',num2str(yr(itd)),') ', ...
                var,'_',num2str(mo(itd)),'.mat'])
            load(infile);
            infile_prev = infile;
        end
        
        % VAR exists because it is in "infile"
        AA = squeeze(VAR(ind(itd),:,:));
        
        AAA = interp2(LON,LAT,AA,lon_rho,lat_rho);
        
        if debug
            disp(['writing AAA(1,1) = ',num2str(AAA(1,1))]);
        end
        
        nc_varput(frcname, ncvarname, scalefactor*AAA + scalefactor2, ...
            [iii-1 0 0], [1 M L]);
        
    end % end of time loop
    
    dt_toc = toc;
    disp(['  took ',num2str(round(dt_toc)),' sec']);
    
end % end of variable loop

