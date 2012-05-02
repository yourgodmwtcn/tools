% Plot vertical line at x
%       [] = linex(x,label,color)

function [] = linex(x,label,color)
    
    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    if size(x,2) == 1, x = x'; end
    
    figure(hFig);
    hold on;
    
    if ~exist('label','var'), label = num2str(x'); end
    if ~exist('color','var'), color = 'k'; end
    
    yax = (get(hAxis,'YLim'));
    xtick = get(hAxis,'XTick');
    %xax = (get(hAxis,'XLim'));
    
    %loglog([x; x]',repmat(yax,length(x),1),'k-','LineWidth',1.5);
    
    for i=1:length(x)       
        plot([x(i) x(i)],yax,'--','LineWidth',1.5,'Color',color);
%        set(hAxis,'XTick',sort([xtick x(i)]));
        %if mod(i,2)
        if ~isempty(label)
            text(double(x(i)),double(yax(end)-(length(label))/2),label,'Rotation',0,'VerticalAlignment','Bottom','FontSize',10,'Color',color);%,'FontWeight','Bold');
        end
        %else
        %    text(x(i),yax(mod(i,2)+2),num2str(1./factor./x(i)),'Rotation',90,'VerticalAlignment','Bottom','FontWeight','Bold');
        %end
    end