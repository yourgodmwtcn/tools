% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       mod_movie(fname, varname, tindices, volume,axis, index, commands)
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end] OR [start count] OR [start] - single snapshot
%           volume - cell array with location of selection edges
%                       {'x' x1 x2; 'y' y1 y2; ...} 
%                    x1,x2 etc. can be numbers (index) or string with
%                    actual dimensional values
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice or
%                 - location in proper units (string; finds nearest point)
%                 - 'mid' will find midpoint of 'axis' for you
%           commands - extra commands to be executed after each plot
%                      (passed to animate.m - optional)

% CHANGELOG:
% Changed to use updated roms_var_grid output                   26 Feb 2012
% Now using roms_ncread_params & roms_tindices                  24 Feb 2012
% Bugfix in read_count - skipped one timestep everytime         23 Feb 2012
% Bugfixes - title index when timesteps(1) ~= 1,                22 Feb 2012
%          - change depth to -ve by default and warn user
% Made the animate call generic - woot! & bugfixes              21 Feb 2012
%   - Imposed same colorbar even when multiple strides are used   
%   - Fixed bug with Esc quitting when multiple strides are used  
%   - Now uses km for x,y axes if applicable
%   - Passes stride number (i) and dt to animate.m
% Extended tindices options and fixed major bug in              14 Feb 2012
%      displaying subsets
% Added 'mid' as possible input for index                       07 Feb 2012
% Fixed 'z' plane sections & added filename check               06 Feb 2012
% Fixed display of single timestep                              26 Jan 2012
% Reducd slab to 100 and fixed many bugs                        26 Jan 2012
%   - titles work perfectly now. Added labels.tmax
%   - throws error if given a 2d simulation and asking for x-y plot
% Changed to mod_movie.m + identifies model from file           26 Jan 2012
% Added slab feature                                            26 Jan 2012
% Added extra commands argument to be passed to animate.m
% Uses NaN's to remove zeros at walls to prevent this from mucking contour plots
% Added labels structure
% Original roms_movie version

function [] = mod_movie(fname, varname, tindices, volume, axis, index, commands)

% fname = find_file(fname);
% if isempty(fname) 
%     error('No file found.');
% else
% 	fprintf('\n Using file: %s\n', fname)
% end

%% model specific setup
gcm = 1;
try 
    ncreadatt(fname,'/','MITgcm_version');
catch ME
    gcm = 0;
end

if gcm
    % fix lowercase - uppercase problem
    if length(varname) == 1
        varname = upper(varname); 
    else
        varname(1) = upper(varname(1));
    end
    
    % Set up grid
    [xax,yax,zax,xunits,yunits] = gcm_var_grid(fname,varname);
    time = double(ncread(fname,'T'))./3600/24;

else
    % set up grid    
    [xax,yax,zax,vol] = roms_extract(fname,varname,volume);
    [~,~,~,time,xunits,yunits] = roms_var_grid(fname,varname);
    time = time./3600/24;
end

%% fix input and get needed info

vinfo  = ncinfo(fname,varname);
dim    = length(vinfo.Size);
slab   = 100; % slab for ncread. read 'n' records in at a time - faster response + save some memory?
midflag = 0;  % 1 if script needs to compute the mid level for plot

if ~exist('tindices','var'), tindices = []; end
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

if strcmp(varname,'Eta') || strcmp(varname,'zeta')
    axis = 'z';
    index = 1;
end

if ~exist('commands','var'), commands = ''; end

%% Plot according to options

labels.tmax = time(tindices(1) + dt*floor((tindices(2)-tindices(1))/dt));
labels.revz  = 0;
labels.tunits = 'days';
labels.dt = dt;
labels.t0 = tindices(1)-1;
vartitle = [varname ' (' ncreadatt(fname,varname,'units') ') | '];

figure;

if strcmp(index,'mid'), midflag = 1; end

for i=0:iend-1
    % if reading data in multiple strides, an escape in the first stride should
    % not results in the next stride getting read / animated.
    if i>0 && strcmp(get(gcf,'currentkey'),'escape')
        return;
    end
    
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt,vol);
    
    labels.time = time(read_start(end):dt:(read_start(end)+(read_count(end))*dt -1)); % read_start(end)-1
    labels.stride = i;
            
    switch axis
        case 'x'
            sliceax = xax; 
            plotx = yax;
            ploty = zax;
            axind = 1;
            labels.xax = ['Y (' yunits ')'];
            labels.yax = 'Z (m)';
            
        case 'y'
            sliceax = yax;
            plotx = xax;
            ploty = zax;
            axind = 2;
            labels.xax = ['X (' xunits ')'];
            labels.yax = 'Z (m)';

        case 'z'
            sliceax = zax;
            plotx = xax;
            ploty = yax;
            axind = 3;
            labels.xax = ['X (' xunits ')'];
            labels.yax = ['Y (' yunits ')'];
            
            if ischar(index) && str2double(index) > 0                
                warning('Changed input depth %s m to -%s m', index, index);
                index = num2str(-1 * str2double(index));
            end
            
            if dim == 3 % catch zeta - free surface elevation
                stride = [1 1 dt];
                midflag = 0;
                index  = 1;
                sliceax(index) = 0;
            end

        otherwise
            error('Invalid axis label.');
    end
    
    %% generic animate call
    
    % given location instead of index
    if midflag, index = num2str((sliceax(1)+sliceax(end))/2); end
    if ischar(index), index = find_approx(sliceax,str2double(index),1); end
    
    % fix title string
    labels.title = [vartitle axis ' = ' sprintf('%5.2f', sliceax(index)) ' m | '];

    % read data
    if dim ~= 3
        read_start(axind) = index;
        read_count(axind) = 1;
    end
    dv = double(squeeze(ncread(fname,varname,read_start,read_count,stride)));  
    
    % take care of walls for mitgcm - NEEDS TO BE CHECKED
    if gcm
        s = size(dv);            
        s(3) = size(dv,3); % correct for single timestep snapshots - in that case s is a 2 element row
        if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
        if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
        
        if axind == 3
            if repnan(dv(:,1,:),0)   == zeros([s(1) 1 s(3)]), dv(:,1,:)   = NaN; end;
            if repnan(dv(:,end,:),0) == zeros([s(1) 1 s(3)]), dv(:,end,:) = NaN; end;
        end
    end
    
    s   = size(dv);            
    if s(1) == 1 || s(2) == 1
        close(gcf);
        error('2D simulation?');
    end
    
    if max(plotx(:)) > 1000
        plotx = plotx/1000;
        labels.xax = [labels.xax ' x 10^3'];
    end
    if max(ploty(:)) > 1000
        ploty = ploty/1000;
        labels.yax = [labels.yax ' x 10^3'];
    end
         
    animate(plotx,ploty,dv,labels,commands,3,0.1);
end
