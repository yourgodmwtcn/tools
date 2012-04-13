% Load chindi variables 
%       [misc] = roms_load_misc(fname)

function [misc] = roms_load_misc(fname)
    
    if ~exist('fname','var'), fname = 'ocean_his.nc'; end
    
    vars = {'nl_tnu2','nl_visc2','rdrg','rdrg2','rho0','R0','Tcoef','Scoef','h','f'};
    
    for i=1:length(vars)
        command = ['misc.' vars{i} '= ncread(''' fname ''',''' vars{i} ''');'];
        try
            eval(command);
        catch ME
            command = ['misc.' vars{i} '= NaN;'];
            eval(command);
        end
    end
    