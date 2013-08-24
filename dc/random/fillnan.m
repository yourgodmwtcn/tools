% Replaces val with NaN's in 'in'
%       [out] = fillnan(in,val)

function [out] = fillnan(in,val,recursive)
    in = double(in);
    mask = in == val;
    out = in;
    out(mask) = NaN;
    
    if nargin == 2
        recursive = 0;
    end
    
    if val == Inf || val == -Inf
        if recursive == 1, return; end
        fprintf('\n Running fillnan recursively for +-Inf \n');
        out = fillnan(out,-1*val,1);
    end
        