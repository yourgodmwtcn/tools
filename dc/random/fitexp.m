% Fits exponential of form fun = @(a,x) a(1) + a(2)*exp(a(3)*x);
% Initial guess a = [0;0;0]. If it doesn't work supply a0.
% Returns a
%           [a] = fitexp(xdata,ydata,a0)

function [a] = fitexp(xdata,ydata,a0)

    if ~exist('a0','var'), a0 = [0;0;0]; end
    fun = @(a1,xdata) a1(1) + a1(2)*exp(a1(3)*xdata);
    [a,~,~,exitflag] = lsqcurvefit(fun,a0,xdata,ydata);
    exitflag
    if exitflag <= 0, error('No fit possible. Change a0.'); end