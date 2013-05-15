% [params] = read_params_from_ini(fname)
% Pass either *.in file or directory

function [params] = read_params_from_ini(fname)

    if isdir(fname)
        if isempty(strfind(fname,'config')), fname = [fname '/config/']; end
        in = ls([fname '/*.in']);
        [~,p] = grep('ININAME == ',[fname in]); 
        % line in p.match must be processed to extract *.nc name
        A = sscanf(char(p.match),' ININAME == %s');
        fname = [fname '/' A];
    end
    info = ncinfo(fname);
    n = length(info.Variables);
    
    for i=1:n
        name = info.Variables(i).Name;
        if strfind(name,'.')
            eval(['params.' name ' = ncread(''' fname ''',''' name ''');']);
        end
    end