% Movie of field 'varname' from 'fname' sliced along 'index' of 'axis'.
%       mod_movie(model,fname, varname, tindices, axis, index, commands)
%           fname - filename
%           varname - variable name
%           tindices - time indices to animate [start (stride) end]
%           axis - axis of slice ('x','y' or 'z')
%           index - index along 'axis' to slice or
%                 - location in proper units (finds nearest point)
%                 - 'mid' will find midpoint of 'axis' for you
%           commands - extra commands to be executed after each plot
%                      (passed to animate.m - optional)

% CHANGELOG:
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

fname = find_file(fname);
if isempty(fname)
    error('No file found.');
else
	fprintf('\n Using file: %s\n', fname)
end

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

if ~exist('tindices','var') || isempty(tindices)
    tindices = [1 Inf];
    dt = 1;
else    
    if length(tindices) == 3
        dt = tindices(2);
        tindices(2) = tindices(3);
        tindices(3) = NaN;
    else
        if length(tindices) == 2
            dt = 1;
        end    
    end
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
vartitle = [varname ' (' ncreadatt(fname,varname,'units') ') | '];

figure;

for i=0:iend-1
    % start and count arrays for ncread : corrected to account for stride
    read_start = ones(1,dim);
    read_end   = Inf(1,dim);
    read_start(end) = tindices(1)+jump*i;
    if i == (iend-1)
        read_end(end) = ceil((tindices(2)-jump*(i))/dt);
    else
        read_end(end) = ceil(jump/dt)-1;%ceil(jump*(i+1)/dt);
    end
    
    labels.time = time(read_start(end):dt:(read_start(end)+read_end(end)*dt-1));
    
    if strcmp(index,'mid'), midflag = 1; end
            
    switch axis
        case 'x'
            % given location instead of index
            if midflag, index = num2str((xax(1)+xax(end))/2); end
            
            if ischar(index), index = find_approx(xax,str2double(index),1); end
            % fix title string
            labels.title = [vartitle axis ' = ' sprintf('%5.2f', xax(index)) ' m | '];
            
            read_start(1) = index;
            read_end(1)   = 1;
            dv = double(squeeze(ncread(fname,varname,read_start,read_end,stride)));            
            %dv  = double(squeeze(var));
            
            % take care of walls
            s   = size(dv);            
            s(3) = size(dv,3); % correct for single timestep snapshots
            if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
            if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
            
            labels.xax = ['Y (' yunits ')'];
            labels.yax = 'Z (m)';
            animate(yax,zax,dv,labels,commands);

        case 'y'
            % given location instead of index
            if midflag, index = num2str((yax(1)+yax(end))/2); end
            if ischar(index), index = find_approx(yax,str2double(index),1); end
            % fix title string
            labels.title = [vartitle axis ' = ' sprintf('%5.2f', yax(index)) ' m | '];
            
            read_start(2) = index;
            read_end(2)   = 1;
            dv = double(squeeze(ncread(fname,varname,read_start,read_end,stride)));
            %dv  = double(squeeze(var));
            
            % take care of walls
            s   = size(dv);
            s(3) = size(dv,3); % correct for single timestep snapshots
            if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:) = NaN; end;
            if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
            
            labels.xax = ['X (' xunits ')'];
            labels.yax = 'Z (m)';
            animate(xax,zax,dv,labels,commands);

        case 'z'
            % given location instead of index            
            if midflag, index = num2str((zax(1)+zax(end))/2); end
            if ischar(index), index = find_approx(zax,str2double(index),1); end
                        
            if dim == 3 % catch zeta
                % dont need to change read_start & read_end
                stride = [1 1 dt];
                % fix title string
                labels.title = [vartitle axis ' = 0 m | '];
            else
                read_start(3) = index;
                read_end(3) = 1;
                stride = [1 1 1 dt];
                % fix title string
                labels.title = [vartitle axis ' = ' sprintf('%5.2f', zax(index)) ' m | '];
            end
            var = ncread(fname,varname,read_start,read_end,stride);
            s1 = size(var);
            dv  = double(squeeze(var));
            clear var;
            
            % take care of walls
            s   = size(dv);            
            s(3) = size(dv,3); % correct for single timestep snapshots
            if s1(1) == 1 || s1(2) == 1
                close(gcf);
                error('2D simulation?');
             end
            if repnan(dv(1,:,:),0)   == zeros([1 s(2) s(3)]), dv(1,:,:)   = NaN; end;
            if repnan(dv(end,:,:),0) == zeros([1 s(2) s(3)]), dv(end,:,:) = NaN; end;
            if repnan(dv(:,1,:),0)   == zeros([s(1) 1 s(3)]), dv(:,1,:)   = NaN; end;
            if repnan(dv(:,end,:),0) == zeros([s(1) 1 s(3)]), dv(:,end,:) = NaN; end;
            
            labels.xax = ['X (' xunits ')'];
            labels.yax = ['Y (' yunits ')'];
            animate(xax,yax,dv,labels,commands);

        otherwise
            error('Invalid axis label.');
    end
end
