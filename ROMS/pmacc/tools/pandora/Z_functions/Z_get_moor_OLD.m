function [M] = Z_get_moor(run,mlon,mlat)
% Z_get_moor.m  3/28/2012  Parker MacCready
%
% super-duper mooring extractor
%
% uses the intrinsic netcdf commands like netcdf.getVar instead of
% SNCTOOLS, and no longer loops through depth.
%
% also the code now automatically extracts ALL 2D and 3D variables at a
% given location.  Note that you can look at the extracted values of
% lon_rho, etc. to check that the interpolation is done correctly.

tic
NT = length(run.his.ncn);
for ttt = 1:NT
    tt = run.his.ncn(ttt);
    tt_str = ['0000',num2str(tt)];
    tt_str = tt_str(end-3:end);
    % get the correct infile
    [infile] = [run.his.dirname,run.his.basename,tt_str,'.nc'];
    if mod(ttt,10)==0; disp(['  infile = ',infile]); end;
    
    ncid = netcdf.open(infile,'NOWRITE');
    
    % get the full listing of variables and their info
    varid_list = netcdf.inqVarIDs(ncid);
    for ii = 1:length(varid_list)
        varid = varid_list(ii);
        [varname_list{ii},xtype,dimids_list{ii},natts] = ...
            netcdf.inqVar(ncid,varid);
    end
    
    % make the indices and weighting matrices for the three grids
    var2get_list = {'lon_rho','lat_rho','lon_u','lat_u', ...
        'lon_v','lat_v','lon_psi','lat_psi'};
    for vv = 1:length(var2get_list)
        var2get = var2get_list{vv};
        tf = strcmp(var2get,varname_list);
        if sum(tf)==1
            % do this if the variable exists
            index = find(tf);
            varname = varname_list{index};
            varid = varid_list(index);
            var = netcdf.getVar(ncid,varid,'double');
            % NOTE that using netcdf.getVar the fields are packed (x,y,z,t) so
            % the lat and lon matrices are packed (x,y), but this is confusing
            % so we will convert EVERYTHING back to (y,x)
            var = var'; % now stored as var(y,x)
            eval([varname,' = var;']);
        else
            % do this if the variable is absent
            disp(['** ',var2get,' not found!']);
        end
    end
    
    % find weighting matrix for the four points closest to the mooring
    % on the rho-, u-, v- and psi-grids
    %
    for gg = 1:4
        switch gg
            case 1
                lon = lon_rho; lat = lat_rho;
            case 2
                lon = lon_u; lat = lat_u;
            case 3
                lon = lon_v; lat = lat_v;
            case 4
                lon = lon_psi; lat = lat_psi;
        end
        ix = find(lon(1,:)<=mlon,1,'last');
        jy = find(lat(:,1)<=mlat,1,'last');
        delx = mlon - lon(1,ix);
        Delx = lon(1,ix+1) - lon(1,ix);
        dely = mlat - lat(jy,1);
        Dely = lat(jy+1,1) - lat(jy,1);
        Rx = delx/Delx; Ry = dely/Dely;
        weight = [(1-Ry)*(1-Rx) (1-Ry)*Rx; Ry*(1-Rx) Ry*Rx];
        switch gg
            case 1
                weightr = weight; ixr = ix; jyr = jy;
            case 2
                weightu = weight; ixu = ix; jyu = jy;
            case 3
                weightv = weight; ixv = ix; jyv = jy;
            case 4
                weightp = weight; ixp = ix; jyp = jy;
        end
        clear ix jy delx Delx dely Dely Rx Ry weight
    end
    clear lon_rho lat_rho lon_u lat_u lon_v lat_v lon_psi lat_psi gg lon lat
    
    % now let's get all the desired variable
    var2get_list = varname_list;
    for vv = 1:length(var2get_list)
        var2get = var2get_list{vv};
        tf = strcmp(var2get,varname_list);
        if sum(tf)==1
            index = find(tf);
            varname = varname_list{index};
            varid = varid_list(index);
            dimids = dimids_list{index};
            if ~isempty(dimids)
                [xname,xlen] = netcdf.inqDim(ncid,dimids(1));
                go_ahead = 0;
                if strcmp(xname(1:2),'xi')
                    go_ahead = 1; % it is an x,y variable
                    switch xname
                        case 'xi_rho'
                            ix=ixr; jy=jyr; weight=weightr;
                        case 'xi_u'
                            ix=ixu; jy=jyu; weight=weightu;
                        case 'xi_v'
                            ix=ixv; jy=jyv; weight=weightv;
                        case 'xi_psi'
                            ix=ixp; jy=jyp; weight=weightp;
                    end
                end
                if go_ahead
                    switch length(dimids)
                        case 2 % XY variables
                            if ttt == 1
                                [yname,ylen] = netcdf.inqDim(ncid,dimids(2));
                                var0 = netcdf.getVar(ncid,varid, ...
                                    [ix-1 jy-1],[2 2],'double');
                                var0 = permute(var0,[2 1]);
                                var = squeeze(sum(sum(weight.*var0,2),1));
                                eval(['M.',varname,' = var;']);
                            end
                        case 3 % XYT variables
                            [yname,ylen] = netcdf.inqDim(ncid,dimids(2));
                            [tname,tlen] = netcdf.inqDim(ncid,dimids(3));
                            var0 = netcdf.getVar(ncid,varid, ...
                                [ix-1 jy-1  0],[2 2 1],'double');
                            var0 = permute(var0,[2 1]);
                            var = squeeze(sum(sum(weight.*var0,2),1));
                            eval(['M.',varname,'(ttt) = var;']);
                        case 4 % XYZT variables
                            [yname,ylen] = netcdf.inqDim(ncid,dimids(2));
                            [zname,zlen] = netcdf.inqDim(ncid,dimids(3));
                            [tname,tlen] = netcdf.inqDim(ncid,dimids(3));
                            var0 = netcdf.getVar(ncid,varid, ...
                                [ix-1 jy-1 0  0],[2 2 zlen 1],'double');
                            var0 = permute(var0,[2 1 3]);
                            var = squeeze(sum(sum(repmat(weight,[1 1 zlen]).*var0,2),1));
                            eval(['M.',varname,'(:,ttt) = var;']);
                    end
                end
            end
        else
            disp(['** ',var2get,' not found!'])
        end
    end
    
    netcdf.close(ncid);
    
end % end of ttt loop
dt = toc;
disp([num2str(round(dt)),' sec for ',num2str(NT),' saves'])

% calculate z-position
M.z_rho = roms_z(M.h*ones(1,NT)',M.zeta',run.grid.cs);
M.z_w = roms_z(M.h*ones(1,NT)',M.zeta',run.grid.csw);
%
% and add a few things for backward compatibility
M.p_lon = mlon; M.p_lat = mlat; M.hh = M.h;

