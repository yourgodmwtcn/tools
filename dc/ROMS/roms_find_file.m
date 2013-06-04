% finds files for type = ini / bry / grd / flt / his / avg
%   [fname] = roms_find_file(dir,type)
% can return a list for last two types. fname does NOT contain input
% directory

function [fname] = roms_find_file(dir,type)
    
    if isempty(strfind(dir,'config')), fname = [dir '/config/']; end
    in = ls([fname '/*.in']);
    
    % files from *.in 
    if strcmpi(type,'ini') || strcmpi(type,'bry') ||  strcmpi(type,'grd')
        fname = ['/config/' grep_in([fname in],type)];
    end
    
    % floats
    if strcmpi(type,'flt')
        fname = [grep_in([fname in],type)];
    end
    
    if strcmpi(type,'his') || strcmpi(type,'avg')
        fname = ls([dir '/*_his*.nc']);
        if isempty(fname)
            fname = ls([dir '/*_avg*.nc']); 
            if strcmpi(type,'his'), disp('Using avg files instead.'); end
        end
        
    end
    

% runs grep on input file
function [str] = grep_in(fname,type)
        [~,p] = grep('ININAME == ',fname); 
        % line in p.match must be processed to extract *.nc name
        str = sscanf(char(p.match),sprintf(' %sNAME == %%s',upper(type)));
        