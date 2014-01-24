function xj = jitter(x,r);

% xj = jitter(x,r);
%
% neil banas, apr 2012

xj = x + r.*randn(size(x));