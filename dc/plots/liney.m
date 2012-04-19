% Plots horizontal line at a given y (can be a vector)
%       [] = liney(y,label,color)

function [] = liney(y,label,color)
    
    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    figure(hFig);
    hold on;
    
    if ~exist('label','var') || isempty(label), label = ' '; end
    if ~exist('color','var'), color = 'k'; end
    
    %yax = (get(hAxis,'YLim'));
    xax = (get(hAxis,'XLim'));
    ytick = get(hAxis,'YTick');    
    
    %loglog([x; x]',repmat(yax,length(x),1),'k-','LineWidth',1.5);
    
    for i=1:length(y)       
        plot(xax,[y(i) y(i)],'--','LineWidth',1.5,'Color',color);
%        set(hAxis,'YTick',sort([ytick y(i)]));
        %if mod(i,2)
            text(double(xax(end)-(length(label)))/2,double(y(i)),label,'Rotation',0,'VerticalAlignment','Bottom','FontSize',10,'Color',color);%,'FontWeight','Bold');
        %else
        %    text(x(i),yax(mod(i,2)+2),num2str(1./factor./x(i)),'Rotation',90,'VerticalAlignment','Bottom','FontWeight','Bold');
        %end
    end