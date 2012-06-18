% Given a cell 'volume' with locations, figures out indices on the axes and
% returns correct axes for the section. 
%           [xax,yax,zax,vol] = roms_extract(fname,varname,volume)
%   volume is cell array of form
%           {'x' 'location1' 'location2';
%            'y' 'location1' 'location2';
%            'z' 'location1' 'location2'};
% locations in strings - in units else provide axis index as number

function [xax,yax,zax,vol] = roms_extract(fname,varname,volume)

    [xax,yax,zax,~,~] = roms_var_grid(fname,varname);
    
    vol = [1 Inf; 1 Inf; 1 Inf]; % default - choose all data
    
    if ~exist('volume','var'), return; end
    if isempty(volume), return; end
    
    % parse volume
    sa = size(volume);
    for i=1:sa(1)
        switch volume{i,1}
            case 'x'
                if ischar(volume{i,2}), volume{i,2} = find_approx(xax,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(xax,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = length(xax); end
                
                xax = xax(volume{i,2}:volume{i,3});
                
                vol(1,1) = volume{i,2};
                vol(1,2) = volume{i,3};
                
            case 'y'
                if ischar(volume{i,2}), volume{i,2} = find_approx(yax,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(yax,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = length(yax); end
                
                yax = yax(volume{i,2}:volume{i,3});
                
                vol(2,1) = volume{i,2};
                vol(2,2) = volume{i,3};
                
            case 'z'
                if ischar(volume{i,2}), volume{i,2} = find_approx(zax,str2double(volume{i,2}),1); end
                if ischar(volume{i,3}), volume{i,3} = find_approx(zax,str2double(volume{i,3}),1); end
                
                if isinf(volume{i,3}), volume{i,3} = length(zax); end
                
                if volume{i,3} < volume{i,2}
                    temp = volume{i,3};
                    volume{i,3} = volume{i,2};
                    volume{i,2} = temp;
                end
                
                zax = zax(volume{i,2}:volume{i,3});
                
                vol(3,1) = volume{i,2};
                vol(3,2) = volume{i,3};
                
            otherwise
                error('Unknown option in volume');
        end
    end