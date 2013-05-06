% [params] = read_params_from_ini(fname)

function [params] = read_params_from_ini(fname)

    info = ncinfo(fname);
    n = length(info.Variables);
    
    for i=1:n
        name = info.Variables(i).Name;
        if strfind(name,'.')
            eval(['params.' name ' = ncread(''' fname ''',''' name ''');']);
        end
    end