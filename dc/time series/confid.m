%%confid.m alpha percent confidence limits for a chi-square variate
%%follows Jenkins and Watts, p. 81 and uses Peter Shaw's chisquare
%%program 
%%nu is number of degrees of freedom
%%for 95% confidence set alpha =.05
%%c. wunsch, may 1995
function [lower,upper]=confid(alpha,nu)
%%get upper tail:

upperv=chisquat(1-alpha/2,nu);

%%get lower tail:
lowerv=chisquat(alpha/2,nu);
lower=nu/upperv;
upper=nu/lowerv;
%%should be sigma^2/S^2 confidence bounds where sigma^2 is true variance
%%check value (J&W) is alpha =.05, nu=19, lower bound is .58
%%upper bound is 2.11


