% Draws horizontal or vertical line and adds tick mark
% called by linex and liney
function [handles] = dcline(ax,x,label,color)

    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    if size(x,2) == 1, x = x'; end
    
    figure(hFig);
    hold on;
    
    if ~exist('label','var'), label = num2str(x'); end
    if ~exist('color','var'), color = 'k'; end
    
    if ax == 'x'
        yax  = get(hAxis,'YLim');
        tickstr = 'XTick';
    else
        xax  = get(hAxis,'XLim');
        tickstr = 'YTick';
    end
    
    tick = get(hAxis,tickstr);
    
    handles = nan(size(x));
    
    for i=1:length(x)
        if ax == 'x'
            handles(i) = plot([x(i) x(i)],yax,'--','LineWidth',2,'Color',color);
            if ~isempty(label)
                text(double(x(i)),double(yax(end)-(length(label))/2),label, ...
                    'Rotation',0,'VerticalAlignment','Bottom','FontSize',10,'Color',color);%,'FontWeight','Bold');
            end
        else
            handles(i) = plot(xax,[x(i) x(i)],'--','LineWidth',2,'Color',color);
            if ~isempty(label)
                text(double(xax(end)-(length(label)))/2,double(x(i)),label, ...
                    'Rotation',0,'VerticalAlignment','Bottom','FontSize',10,'Color',color);%,'FontWeight','Bold');
            end
        end
        set(hAxis,tickstr,sort(unique([tick x(i)])));
    end