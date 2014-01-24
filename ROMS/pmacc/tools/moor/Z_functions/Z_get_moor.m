function [M] = Z_get_moor(run,mlon,mlat)
% Z_get_moor.m  8/19/2013  Parker MacCready
%
% super-duper mooring extractor
%
% uses the intrinsic netcdf commands like netcdf.getVar instead of
% SNCTOOLS, and no longer loops through depth.
%
% also the code now automatically extracts ALL 2D and 3D variables at a
% given location.  Note that you can look at the extracted values of
% lon_rho, etc. to check that the interpolation is done correctly.
%
% 6/6/2012 sng included extraction of ocean_time, preallocated some
% variables, and removed extraneous variables
%
% 1/14/2013 Now accepts vectors of mlon, mlat locations, and returns a
% struct array M of the same size.
%
% 8/19/2013 Replaced call to roms_z with call to Z_s2z (ae end).

% get the year of the run
yeari = datestr(run.his.nctime(1),'yyyy');
yeari = str2double(yeari);

% get info on what variables are available from the first file
ttt = 1;
tt = run.his.ncn(ttt);
tt_str = ['0000',num2str(tt)];
tt_str = tt_str(end-3:end);
% get the correct infile
[infile] = [run.his.dirname,run.his.basename,tt_str,'.nc'];

ncid = netcdf.open(infile,'NOWRITE');

% get the full listing of variables and their info
varid_list = netcdf.inqVarIDs(ncid);
varname_list = cell(1,length(varid_list));
dimids_list = cell(1,length(varid_list));
for ii = 1:length(varid_list)
    varid = varid_list(ii);
    [varname_list{ii},~,dimids_list{ii},~] = ...
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

netcdf.close(ncid);

% find weighting matrix for the four points closest to the mooring
% on the rho-, u-, v- and psi-grids
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
    
    for mmm = 1:length(mlon)
        
        ix = find(lon(1,:)<=mlon(mmm),1,'last');
        jy = find(lat(:,1)<=mlat(mmm),1,'last');
        delx = mlon(mmm) - lon(1,ix);
        Delx = lon(1,ix+1) - lon(1,ix);
        dely = mlat(mmm) - lat(jy,1);
        Dely = lat(jy+1,1) - lat(jy,1);
        Rx = delx/Delx; Ry = dely/Dely;
        weight = [(1-Ry)*(1-Rx) (1-Ry)*Rx; Ry*(1-Rx) Ry*Rx];
        switch gg
            case 1
                weightr(mmm,:,:) = weight; ixr(mmm) = ix; jyr(mmm) = jy;
            case 2
                weightu(mmm,:,:) = weight; ixu(mmm) = ix; jyu(mmm) = jy;
            case 3
                weightv(mmm,:,:) = weight; ixv(mmm) = ix; jyv(mmm) = jy;
            case 4
                weightp(mmm,:,:) = weight; ixp(mmm) = ix; jyp(mmm) = jy;
        end
        clear ix jy delx Delx dely Dely Rx Ry weight
        
    end % end of mmm loop
    
end % end of the gg loop
clear lon_rho lat_rho lon_u lat_u lon_v lat_v lon_psi lat_psi gg lon lat

%loop through the input files
tic
NT = length(run.his.ncn);
for ttt = 1:NT
    tt = run.his.ncn(ttt);
    tt_str = ['0000',num2str(tt)];
    tt_str = tt_str(end-3:end);
    % get the correct infile
    [infile] = [run.his.dirname,run.his.basename,tt_str,'.nc'];
    if mod(ttt,100)==0; disp(['  infile = ',infile]); end;
    
    ncid = netcdf.open(infile,'NOWRITE');
    
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
                [xname,~] = netcdf.inqDim(ncid,dimids(1));
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
                                for mmm = 1:length(mlon)
                                    var0 = netcdf.getVar(ncid,varid, ...
                                        [ix(mmm)-1 jy(mmm)-1],[2 2],'double');
                                    var0 = permute(var0,[2 1]);
                                    this_weight = squeeze(weight(mmm,:,:));
                                    var = squeeze(sum(sum(this_weight.*var0,2),1));
                                    eval(['M(mmm).',varname,' = var;']);
                                end % end of mmm loop
                            end
                        case 3 % XYT variables
                            for mmm = 1:length(mlon)
                                var0 = netcdf.getVar(ncid,varid, ...
                                    [ix(mmm)-1 jy(mmm)-1  0],[2 2 1],'double');
                                var0 = permute(var0,[2 1]);
                                this_weight = squeeze(weight(mmm,:,:));
                                var = squeeze(sum(sum(this_weight.*var0,2),1));
                                eval(['M(mmm).',varname,'(ttt) = var;']);
                            end % end of mmm loop
                        case 4 % XYZT variables
                            for mmm = 1:length(mlon)
                                [~,zlen] = netcdf.inqDim(ncid,dimids(3));
                                var0 = netcdf.getVar(ncid,varid, ...
                                    [ix(mmm)-1 jy(mmm)-1 0  0],[2 2 zlen 1],'double');
                                var0 = permute(var0,[2 1 3]);
                                this_weight = squeeze(weight(mmm,:,:));
                                var = squeeze(sum(sum(repmat(this_weight,[1 1 zlen]).*var0,2),1));
                                eval(['M(mmm).',varname,'(:,ttt) = var;']);
                            end % end of mmm loop
                    end
                end
            end
            %get the ocean time as a cross check to run.his.nctime
            if strcmp(varname,'ocean_time');
                ot = netcdf.getVar(ncid,varid,'double');
                %convert the ocean time to matlab date
                for mmm = 1:length(mlon)
                    M(mmm).td(ttt) = (ot./(24*60*60))+datenum([yeari,1,1]);
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
for mmm = 1:length(mlon)
    if 0
        M(mmm).z_rho = roms_z(M(mmm).h*ones(1,NT)',M(mmm).zeta',run.grid.cs);
        M(mmm).z_w = roms_z(M(mmm).h*ones(1,NT)',M(mmm).zeta',run.grid.csw);
    else
        % NOTE: 8/19/2013 need to replace call to roms_z with call to Z_s2z
        S = Z_get_S(infile);
        [M(mmm).z_rho,M(mmm).z_w] = ...
            Z_s2z(M(mmm).h*ones(1,NT)',M(mmm).zeta',S);
    end
end

