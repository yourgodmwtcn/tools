% Plots horizontal line at a given y (can be a vector)
%       [handles] = liney(y,label,color)

function [handles] = liney(y,label,color)
    
    if ~exist('label','var'), label = []; end
    if ~exist('color','var'), color = 'k'; end
    handles = dcline('y',y,label,color);