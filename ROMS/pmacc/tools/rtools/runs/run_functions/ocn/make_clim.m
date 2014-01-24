function make_clim(indir0,outdir,gridfile,td_to_get,ts_to_get,S,clmname)
% make_clim.m  12/18/2012  Parker MacCready
%
% This creates a NetCDF file to be used for ROMS forcing.  It is designed
% to read in the files created by:
% tools_data/mossea_forcing_data/ocn/NCOM_global/preprocess.m
%
% It makes: a ROMS input climatology file
% with information on T, S, u, v, ubar, vbar and zeta
% sng 12/27/2011 fixed ssh such that missing values have 0 rather than NaN
% PM 5/24/2012 recoded to have cleaner handling across year boundaries

% get vectors of information describing the data locations
[INV] = inventory(indir0,td_to_get);
nts = length(ts_to_get);

% make the output directory if needed
if ~exist(outdir,'dir'); mkdir(outdir); end;

outfile = [outdir,clmname];
clim_netcdf_new(gridfile,outfile,S.N);

h = nc_varget(gridfile,'h');
lon_rho = nc_varget(gridfile,'lon_rho');
lat_rho = nc_varget(gridfile,'lat_rho');

varname_list={'s3d';'t3d';'u3d';'v3d'};
for ivv = 1:length(varname_list) % start of VARNAME LOOP
    tic
    varname = varname_list{ivv};
    [V,AAlon,AAlat,AAmask] = ...
        clim_netcdf(indir0,varname,gridfile);
    [M,L] = size(AAlon);
    hh = interp2(lon_rho,lat_rho,h,AAlon,AAlat);
    disp(['    - writing ',V.invarname,' to ',V.ncvarname]);
    
    clim_netcdf_time(outfile,V,nts);
    if V.do_addtimevar
        nc_varput(outfile, V.nctimename, ts_to_get);
    end
    
    switch varname
        case 'u3d'
            ubar_mat = NaN * ones(nts,M,L);
        case 'v3d';
            vbar_mat = NaN * ones(nts,M,L);
    end
    
    varstruct.Name = V.ncvarname;
    varstruct.Dimension = {V.nctimename, 's_rho', ...
        V.ncetaname, V.ncxiname};
    long_name = V.nclongname;
    units = V.ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(outfile, varstruct);
    fn = [indir0,INV.indir(1,:),'/',varname,'.nc'];
    indepth = nc_varget(fn, 'depth');
    AN = length(indepth);
    
    for tt = 1:length(INV.td) % start of TIME LOOP
              fn = [indir0,INV.indir(tt,:),'/',varname,'.nc'];
  
        if tt == 1
            % NCOM grid
            Alon = nc_varget(fn, 'lon');
            Alat = nc_varget(fn, 'lat');
            [ALON,ALAT] = meshgrid(Alon,Alat);
            % ROMS grid
            [z_rho,z_w] = Z_s2z(hh,0*hh,S);
            DZ = diff(z_w,1,1);
        end
                
        % reinitialize the storage arrays
        % AA is the empty matrix on the ROMS grid in which to put
        % the field from this time step, BUT it has the vertical grid
        % of the NCOM file
        AA = NaN * ones(AN,M,L);
        % AAA is the empty matrix on the ROMS grid in which to put
        % the field from this time step, and it has the vertical grid
        % of the ROMS model.
        AAA = NaN * ones(S.N,M,L);
        
        % load the variable field at this time step
        clear A
        % INTERPOLATE HORIZONTALLY TO ROMS GRIDS (fills whole domain,
        % if there is ANY data on that level)
        A = nc_varget(fn, V.invarname,[INV.index(tt)-1 0 0 0],...
            [1 -1 -1 -1]);
        for zz = 1:AN % start of Z-LOOP
            this_A = squeeze(A(zz,:,:));
            AA(zz,:,:) = interp2(ALON,ALAT,this_A,AAlon,AAlat);
        end
        
        % INTERPOLATE VERTICALLY TO ROMS S-SURFACES
        fillv = 0;
        this_z_ncom = -indepth;
        for jj = 1:M % start of XY-LOOP
            for ii = 1:L
                % only do interpolation on non-masked points
                if AAmask(jj,ii)
                    this_AA = AA(:,jj,ii);
                    this_z_roms = z_rho(:,jj,ii);
                    this_AAA = interp1q(this_z_ncom, this_AA, ...
                        this_z_roms);
                    AAA(:,jj,ii) = this_AAA;
                else
                    AAA(:,jj,ii) = fillv;
                end 
            end 
        end 
        
        switch varname
            case 'u3d'
                ubar_mat(tt,:,:) = squeeze(sum(AAA.*DZ,1)) ./ hh;
                if tt == nts
                    tempufile = [outdir 'temporary_ubar_storage.mat'];
                    save(tempufile,'ubar_mat')
                end
            case 'v3d'
                vbar_mat(tt,:,:) = squeeze(sum(AAA.*DZ,1)) ./ hh;
                if tt == nts
                    tempvfile = [outdir 'temporary_vbar_storage.mat'];
                    save(tempvfile,'vbar_mat')
                end
        end
        
        % write the field at one time level
        [nz,nx,ny]=size(AAA);
        nc_varput(outfile, V.ncvarname, AAA, [tt-1 0 0 0], [1 nz nx ny]);
        clear nz nx ny AAA
        
    end % end of TIME LOOP
    delta_t = toc;
    disp(['      ',num2str(round(delta_t)),' sec for ', ...
        num2str(nts),' times'])
end

% ************************ SSH **********************************
varname = 'ssh';
[V,AAlon,AAlat,AAmask] = ...
    clim_netcdf(indir0,varname,gridfile);
[M,L] = size(AAlon);
disp(['    - writing ',V.invarname,' to ',V.ncvarname]);
    
clim_netcdf_time(outfile,V,nts);
nc_varput(outfile, V.nctimename, ts_to_get);

varstruct.Name = V.ncvarname;
varstruct.Dimension = {V.nctimename, V.ncetaname, V.ncxiname};
long_name = V.nclongname;
units = V.ncunits;
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(outfile, varstruct);

for tt = 1:length(INV.td) % start of TIME LOOP
        fn = [indir0,INV.indir(tt,:),'/',varname,'.nc'];

    if tt == 1
        % these are the NCOM lat and lon for this variable
        Alon = nc_varget(fn, 'lon');
        Alat = nc_varget(fn, 'lat');
        [ALON,ALAT] = meshgrid(Alon,Alat);
    end
        
    % reinitialize the storage arrays
    AA = NaN * ones(M,L); % ROMS grid
    % load the variable field at this time step
    clear A
    % INTERPOLATE HORIZONTALLY TO ROMS GRIDS (fills whole domain)
    A = nc_varget(fn, V.invarname,[INV.index(tt)-1 0 0],[1 -1 -1]);
    AA = interp2(ALON,ALAT,A,AAlon,AAlat);
    %fill any missing values remaining with 0
    AA(isnan(AA)) = 0;
    % write the field at one time level
    [nx,ny]=size(AA);
    nc_varput(outfile, V.ncvarname, AA, [tt-1 0 0], [1 nx ny]);
    clear nx ny AA
    
end % end of TIME LOOP

% ************************ END SSH ******************************

varname_list={'ubar';'vbar'};
for ivv = 1:length(varname_list) % start of VARNAME LOOP
    varname = varname_list{ivv};
    [V,AAlon,AAlat,AAmask] = ...
        clim_netcdf(indir0,varname,gridfile);
    [M,L] = size(AAlon);
    disp(['    - writing ',V.invarname,' to ',V.ncvarname]);
    
    clim_netcdf_time(outfile,V,nts);
    if V.do_addtimevar
        nc_varput(outfile, V.nctimename, ts_to_get);
    end
    
    varstruct.Name = V.ncvarname;
    varstruct.Dimension = {V.nctimename, V.ncetaname, V.ncxiname};
    long_name = V.nclongname;
    units = V.ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(outfile, varstruct);
    
    switch varname
        case 'ubar'
            load(tempufile)
            nc_varput(outfile,V.ncvarname,ubar_mat);
            delete(tempufile);
        case 'vbar'
            load(tempvfile)
            nc_varput(outfile,V.ncvarname,vbar_mat);
            delete(tempvfile);
    end
    
end

disp('  * DONE writing climatology file')
