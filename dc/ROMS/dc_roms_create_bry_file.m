% Function that takes in structure, writes out appropriate grid
% information to boundary conditions file.
% inspired by d_obc_roms2roms.m

function [] = dc_roms_create_bry_file(S)

	OBC.west  = S.boundary(1);
    OBC.east  = S.boundary(2);
    OBC.south = S.boundary(3);
    OBC.north = S.boundary(4);
    
    VarGrd = {'spherical',                                                ...
          'Vtransform', 'Vstretching',                                ...
          'theta_s', 'theta_b', 'Tcline', 'hc',                       ...
          's_rho', 'Cs_r', 's_w', 'Cs_w'};
    
    % write out to file
    for ii=1:length(VarGrd)
        if isfield(S,VarGrd{ii})
            nc_write(S.ncname,VarGrd{ii},S.(VarGrd{ii}));
        end
    end
    
    % write out grid for each boundary
    bry = {'west','east','south','north'};
    out_var = {'u','v','rho'};
    
    if ~S.spherical
        ax = {'x','y'};
    else
        ax = {'lon','lat'};
    end
    
    for ii=1:length(bry)
        % if needed
        if OBC.(bry{ii})
            for jj=1:length(out_var)
                for kk=1:length(ax)
                    part_name = [ax{kk} '_' out_var{jj}];
                    full_name = [ax{kk} '_' out_var{jj} '_' bry{ii}];

                    switch(ii)
                        case 1 % west
                            var = S.(part_name)(1,:);
                        case 2 % east
                            var = S.(part_name)(end,:);
                        case 3 % south
                            var = S.(part_name)(:,1);
                        case 4 % north  
                            var = S.(part_name)(:,end);
                    end               

                    if size(var,1) == 1, var = var'; end

                    nc_write(S.ncname,full_name,var);     
                end
            end
        end                
    end
        
    