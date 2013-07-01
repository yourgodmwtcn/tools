% returns indices where a(vector) = b (vector). Biased high
%       [out] = vecfind(a,b)

function [out] = vecfind(a,b)
    out = arrayfun(@(x) find(a > x,1,'first'), b );