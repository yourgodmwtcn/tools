% Makes animation (default using contourf). Assumes last dimension is to be looped by default. 
% Else specify index. Allows browsing.
%       [mm_instance,handles] = animate(xax,yax,data,labels,commands,index,pausetime)
%           xax,yax - x,y axes (both optional) - can be empty []
%           data - data to be animated - script squeezes data
%
%           labels - Structure with following fields for labeling plot (optional - can be [])
%               > title
%               > xax
%               > yax
%               > revz - reverse yDir?
%               > time - where to pull time info for titles
%               > tunits - units for time vector labels.time
%               > tmax  - maximum value of time vector : needed when only
%                         part of data is being plotted at a time. set by
%                         default to max(labels.time) if not defined
%               > dar - correct display aspect ratio when displaying lat/lon plots (using Parker's code)
%               > stride - which stride are we on when called by mod_movie? (default: 0)
%               > dt - time step interval from mod_movie (default: 1)
%               > t0 - inital timestep from mod_movie
%
%           commands - custom commands to execute after plot or one of below. (in string, optional)
%                    - separate commands using ; 
%                    - Built-in options are: 
%                          > nocaxis - non-constant colorbar
%                          > pause   - start paused
%                          > lab_cmap - LAB space colormap
%                          > pcolor / contour / contourf - imagescnan is default
%                          > movieman - make movie using Ryan's movieman code
%                          > topresent - tweaks image to make it better for saving (bigger fonts, reduced axis tick marks etc.)
%
%           index - dimension to loop through (optional)
%           pausetime - pause(pausetime) (optional)
%
% USAGE:
%       animate(data)
%       animate(data,commands)
%       animate(xax,yax,data)
%       animate([],[],data,...
%       animate(xax,yax,data,commands)
%               eg: animate(xax,yax,data,'pcolor;nocaxis')
%
% BROWSE:
%       - first space pauses, second space resumes, remaining spaces *play*
%       - arrowkeys *pause* and navigate always
%       - Esc to quit

function [mm_instance,handles] = animate(xax,yax,data,labels,commands,index,pausetime)

    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    %% figure out input
    narg = nargin;
    
    if strcmp(get(gcf,'currentkey'),'escape')
        warning('Previous ESC detected. Opening new figure.');
        figure;
    end
        
    switch narg
        case 1,
            data = squeeze(xax);
            s    = size(data);
            xax  = 1:s(1); 
            yax  = 1:s(2);
            
        case 2,
            data = xax;
            commands = yax;
            xax = [];
            yax = [];            
            
        case 4,
            if ischar(class(labels))
               commands = labels;
               labels = [];
            end
    end
    
    data = squeeze(data);
    s = size(data);
    if length(s) == 2, s(3) = 1; end
       
    if narg <= 5 || ~exist('index','var')
        stop = size(data,3);
        index = length(s); % set last dimension to loop by default
    else
        stop = s(index);
    end
    
    if narg ~= 7 || ~exist('pausetime','var')
        pausetime = 0.2;
    end
    
    if ~exist('labels','var') || isempty(labels)
        labels.title = '';
        labels.xax   = '';
        labels.yax   = '';
        labels.revz  = 0;
        labels.time  = [];
        labels.tmax  = size(data,3);
        labels.stride = 0;
        labels.dt = 1;     
        labels.dar = 0;
    end
    
    if ~isfield(labels,'tmax'), labels.tmax = labels.time(end); end
    if ~isfield(labels,'dt'), labels.dt = 1; end
    if ~isfield(labels,'stride'), labels.stride = 0; end
    if ~isfield(labels,'t0'), labels.t0 = 0; end
    
    if ~exist('commands','var'), commands = '';     end
    
    if isempty(xax), xax = [1:s(1)]'; end;
    if isempty(yax), yax = [1:s(2)]'; end;
    
    % make 2D axis variables if they aren't already
    if isrow(xax), xax = xax'; end
    if isrow(yax), yax = yax'; end
    
    if isvector(xax), xax = repmat(xax ,[1 s(2)]); end
    if isvector(yax), yax = repmat(yax',[s(1) 1]); end
    
       
    %% processing
  
    if stop == 1, 
        if ~labels.isDir, warning('Only one time step.'); end
        plotdata = double(squeeze(data)); % shiftdim screws up single snapshots
    else        
        plotdata = double(squeeze(shiftdim(data,index)));
    end    
    
    if labels.stride == 0
        datamax = nanmax(plotdata(:));
        datamin = nanmin(plotdata(:));
    else
        clim = caxis;
    end
    
    hfig = gcf;
    ckey = '';
    button = [];
    pflag = 0;
    spaceplay = 1; % if 1, space pauses. if 0, space plays
    
    %% parse options
    cmds = {'nocaxis','pcolor','contourf','contour','pause','fancy_cmap','movieman','topresent'};
    flags = zeros(1,length(cmds));
    if ~isempty(commands),
        [flags, commands] = parse_commands(cmds,commands);
    end
    
    plotflag = sum([2 3 4] .* flags(2:4));
    if flags(5) && stop ~= 1, spaceplay = 0; fprintf('\n Hit a key to advance frame. \n\n'); end
    
    fancy_map = flipud(cbrewer('div', 'RdYlGn', 32));
    
    if flags(6) % Build colormap
        radius = 50;
        num = 40;
        theta = linspace(0, pi/2, num).';
        a = radius * cos(theta);
        b = radius * sin(theta);
        L = linspace(0, 100, num).';
        Lab = [L, a, b];
        fancy_map = applycform(Lab, makecform('lab2srgb'));
    end
    if flags(7) || flags(8), set(gcf,'Position',[0 0 1600 900]); flags(8)=1; end % maximize for recording + activate topresent

    set(gcf,'Renderer','zbuffer'); % performance!
    i=0;
    while i<=stop-1

        % navigation
        if strcmp(ckey,'space') && isempty(button)
            spaceplay = ~spaceplay;
        end

        pflag = ~spaceplay;
        
        if pflag,
            [~,~,button] = ginput(1);
            figure(gcf);
            if button == 32, spaceplay = 1; end % resumes when paused
            if button == 27, break; end % exit when Esc is pressed.
        else
            pause(0.02);%(pausetime);
        end  
        
        ckey = get(gcf,'currentkey');% end
        
        % navigate : other keys move forward
        if strcmp(ckey,'leftarrow') || strcmp(ckey,'downarrow') || button == 28 || button == 31 || button == 3
            pflag = 1; spaceplay = 0;
            i = i-2;
        else if strcmp(ckey,'rightarrow') || strcmp(ckey,'uparrow') || button == 29 || button == 30 || button == 1
                pflag = 1; spaceplay = 0;
            end
        end
        
        if strcmp(ckey,'escape')
            break
        end
        
        i=i+1;
        if i < 1, i = 1; end
        
        %% Plot
        hold off; % just in case
        switch plotflag
            case 2
                try
                    handles.h_plot = pcolorcen(xax,yax,plotdata(:,:,i));
                catch ME
                    handles.h_plot = pcolor(xax,yax,plotdata(:,:,i));
                end
            case 3
                [~,handles.h_plot] = contourf(xax,yax,plotdata(:,:,i),linspace(datamin,datamax,25)); shading flat
            case 4
                set(gcf,'Renderer','painters');
                clf;
                [C,handles.h_plot] = contour(xax,yax,plotdata(:,:,i),linspace(datamin,datamax,30),'Color','k');
                format short
                clabel(C,handles.h_plot,'FontSize',9);
            otherwise
                try
                    handles.h_plot = imagescnan(xax,yax,plotdata(:,:,i));
                    set(gca,'yDir','normal');
                catch ME
                    handles.h_plot = imagesc(xax,yax,plotdata(:,:,i));
                    set(gca,'yDir','normal');
                end
        end
        
        % square axis if appropriate
        if abs((max(xax)-min(xax)) - (max(yax)-min(yax))) < 1
            axis square;
        end
        
        % colorbar
        if plotflag ~=4 
            if ~flags(1)
                if labels.stride > 0
                    caxis(clim);
                else
                    if datamax ~= datamin,caxis([datamin datamax]); end
                end
            end
        end        
        shading flat;
        colormap(fancy_map);
        if plotflag ~= 4, handles.h_cbar = colorbar;  end
        
        % fix display aspect ratio for lat,lon plots
        if labels.dar, Z_dar; end
        
        % labels
        if labels.revz, revz; end;
        if isempty(labels.time)
            addtitle = [' t instant = ' num2str(labels.t0+i+(labels.dt-1)*(i-1)+100*labels.stride)];
        else
            addtitle = [' t = ' sprintf('%.2f/%.2f ', labels.time(i),labels.tmax) ...
                   labels.tunits];  
            if flags(8) == 0 % don't want t instant in final plots
                addtitle = [addtitle ' (instant = ' num2str(labels.t0+i+(labels.dt-1)*(i-1)+100*labels.stride) ')'];
            end
        end
        
        % center colorbar
        if flags(7) || flags(8)
            [cmin,cmax] = caxis;
            if cmax*cmin < 0 % make colorbar symmetric about zero
                if cmax > abs(cmin)
                    cmin = -abs(cmax);
                else
                    cmax = abs(cmin);
                end
            end
            caxis([cmin cmax]);
        end
        
        eval(commands); % execute custom commands
        beautify;
        % make fonts bigger for presentation / movie plots
        if flags(8), fontSize = 20; else fontSize = 12; end
        
        handles.h_title = title([labels.title addtitle],'FontSize',fontSize);
        xlabel(labels.xax,'FontSize',fontSize);
        ylabel(labels.yax,'FontSize',fontSize);
        set(gca,'FontSize',fontSize);
        
        if flags(7)  
            if isempty(labels.mm_instance)
                labels.mm_instance = mm_setup;
                labels.mm_instance.pixelSize = [1600 900];
                labels.mm_instance.outputFile = 'mm_output.avi';
                labels.mm_instance.ffmpegArgs = '-q:v 1 -g 1';
                labels.mm_instance.InputFrameRate = 3;
                labels.mm_instance.frameRate = 3;
            end
            set(gcf,'Renderer','zbuffer');
            mm_addFrame(labels.mm_instance,gcf);
            mm_instance = labels.mm_instance;             
        else
            mm_instance = [];
        end
        
        jframe = get(gcf,'JavaFrame');
        jframe.setMaximized(true);
        
    end