% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       mod_movie(model,fname, varname, tindices, axis, index, commands)
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end] OR [start count] OR [start] - single snapshot
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice or
%                 - location in proper units (finds nearest point)
%                 - 'mid' will find midpoint of 'axis' for you
%           commands - extra commands to be executed after each plot
%                      (passed to animate.m - optional)

% CHANGELOG:
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
% Reducd jump to 100 and fixed many bugs                        26 Jan 2012
%   - titles work perfectly now. Added labels.tmax
%   - throws error if given a 2d simulation and asking for x-y plot
% Changed to mod_movie.m + identifies model from file           26 Jan 2012
% Added stride feature                                          26 Jan 2012
% Added extra commands argument to be passed to animate.m
% Uses NaN's to remove zeros at walls to prevent this from mucking contour plots
% Added labels structure
% Original roms_movie version

function [] = mod_movie(fname, varname, tindices, axis, index, commands)

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
    
    % set time variable
    time_var = 'T';
else
    % set up grid    
    [xax,yax,zax,xunits,yunits] = roms_var_grid(fname,varname);
    
    % convert matrices to vectors
    xax = xax(:,1);
    yax = yax(1,:);
    if ~isempty(zax), zax = zax(:,1,1); end
    
    % set time variable
    time_var = 'ocean_time';
end

%% fix input and get needed info

vinfo  = ncinfo(fname,varname);
dim    = length(vinfo.Size);
jump   = 100; % jump for ncread. read 'n' records in at a time - faster response + save some memory?
midflag = 0;  % 1 if script needs to compute the mid level for plot

if ~exist('tindices','var') || isempty(tindices)
    tindices = [1 Inf];
    dt = 1;
else
    switch length(tindices)
        case 1
            dt = 1;
            tindices(2) = tindices(1);
            
        case 2
            dt = 1;
        
        case 3
            dt = tindices(2);
            tindices(2) = tindices(3);
            tindices(3) = NaN;
    end
end

if tindices(2) < tindices(1)
    tindices(2) = tindices(1)+tindices(2);
end
if isinf(tindices(2)), tindices(2) = vinfo.Size(end); end

stride = [1 1 1 dt];
    
if (tindices(2)-tindices(1)) == 0
    iend = 1;
    dt = tindices(2);
else
    iend   = ceil((tindices(2)-tindices(1))/jump/dt);
end

if strcmp(varname,'Eta') || strcmp(varname,'zeta')
    axis = 'z';
    index = 1;
end

if ~exist('commands','var'), commands = ''; end

%% Plot according to options
time = double(ncread(fname,time_var))./3600/24;
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
    if i>0 & strcmp(get(gcf,'currentkey'),'escape')
        return;
    end
    
    % start and count arrays for ncread : corrected to account for stride
    read_start = ones(1,dim);
    read_count = Inf(1,dim);
    
    if i == (iend-1)
        read_count(end) = ceil((tindices(2)-jump*(i))/dt);
    else
        read_count(end) = ceil(jump/dt)-1;%ceil(jump*(i+1)/dt);
    end
    
    if i == 0
        read_start(end) = tindices(1);
    else
        read_start(end) = jump*i + 1;
    end
    
    if (iend-1) == 0, read_count(end) = ceil((tindices(2)-tindices(1))/dt)+1; end 
    
    labels.time = time(read_start(end):dt:(read_start(end)+read_count(end)*dt -1)); % read_start(end)-1
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
    
    % take care of walls for mitgcm
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
         
    animate(plotx,ploty,dv,labels,commands);
end
