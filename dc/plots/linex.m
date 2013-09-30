% Plot vertical line at x
%       [handles] = linex(x,label,color)

function [handles] = linex(x,label,color)
    
    if ~exist('label','var'), label = num2str(x'); end
    if ~exist('color','var'), color = 'k'; end
    
    handles = dcline('x',x,label,color);