% adds string to legend
%       [] = addlegend(handle,str)
% from http://www.mathworks.com/support/solutions/en/data/1-181SJ/?solution=1-181SJ

function [] = addlegend(handle,str,loc)

    if ~exist('loc','var') || isempty(loc), loc = 'NorthEast'; end
    legend('show');
    
    % get legend handles
    [~,~,outh,outm] = legend;
    
    % remove default title
    if length(outh) == 1 && strcmp(outm{1},'data1')
        outh = [];
        outm{1} = str;
        legh = legend([outh;handle],outm{:},'');
        return;
    end
    
    legh = legend([outh;handle],outm{:},str);
    set(legh,'Location',loc);
    