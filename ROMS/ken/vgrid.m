%   Experiment with the the ROMS 3 vertical stretching functions
%   Call as
%       vgrid(vxform,vstretch,[param],hc,h,n,plot)
%   where
%       vxform = 1      original formulation
%              = 2      new formulation    (new default)
%       vstretch    = 1 original stretching function
%                   = 2 Shchepetkin (2005) function
%                   = 3 Geyer function
%                   = 4 Shchepetkin(2010) function: New default
%       param = constants for the fit
%               for vstretch = 1:  theta_s, theta_b
%               for vstretch = 2:  theta_s, theta_b, alpha, beta
%               for vstretch = 3:  alpha, beta,  hscale
%               for vstretch = 4:  theta_s, theta_b
%       hc = tcline (m)
%       h = water depth (m)
%       n = # of s points
%       plot = 1 - make plots
%
%       Assumes zeta = 0
%       Theta_s -> surface
%       Theta_b -> bottom
%       hc -> approximate boundary layer depth where lesser spacing is needed

%       K. Brink 1/11/11  8/29/2011

function [s,z] =vgrid(vxform,vstretch,param,hc,h,n,plots)
    
    if ~exist('plots','var'), plots = 0; end

ds = 1/n;
s = ds*(0.5 + (0:(n-1))) -1.;

if round(vstretch) == 1
    ths = param(1);
    thb = param(2);
    cs = (1-thb)*(sinh(s*ths)/sinh(ths)) + thb*(-0.5+0.5*tanh(ths*(s+0.5))/tanh(0.5*ths));
elseif round(vstretch) ==2
    ths = param(1);
    thb = param(2);
    alpha = param(3);
    beta = param(4);
    csur = (1-cosh(ths*s))/(cosh(ths)-1);
    cbot = -1 + sinh(thb*(s+1))/sinh(thb);             %  ????
    cw = ((s+1).^alpha).*(1 + (alpha/beta)*(1 - (s+1).^beta));
    cs = cw.*csur+(1-cw).*cbot;        
elseif round(vstretch) == 3
    ths = param(1);
    thb = param(2);
    hscale = param(3);                  %   use hscale = 3
    csur = -log(cosh(hscale*abs(s).^ths))/log(cosh(hscale));
    cbot = (log(cosh(hscale*(s+1).^thb))/log(cosh(hscale)))-1;
    cw = 0.5*(1-tanh(hscale*(s + 0.5)));
    cs = cw.*cbot + (1-cw).*csur;
    
    
elseif round(vstretch) == 4
    ths = param(1);
    thb = param(2);
    
    %cs = 0*s;              %    Need to fix this!
    if ths > 0
        cs = (1 - cosh(ths*s))/(cosh(ths) - 1);
    else
        cs = -s.*s;
    end
    if thb > 0
        cs = (exp(thb*cs)-1)/(1 - exp(-thb));
    end
    
  
else
    cs = 0*s;
end 
        
        
if round(vxform) == 1
    zz = hc*s +(h-hc)*cs;
    z = zz;
    
elseif round(vxform) == 2
    zz = (hc*s + h*cs)/(hc+h);
    z = h*zz;
else
    z = 0*s;
end

    dz = NaN*z;
    dz(2:end-1) = (z(3:end)-z(1:end-2))/2;
    dzmax = max(dz)
    dzmin = min(dz)

if plots
    figure
    plot(s,z,'+')
    ylabel('z (m)')
    if round(vxform) == 1
        title(['Original form, vstretch = ' int2str(vstretch)])
    else
         title(['New form, vstretch = ' int2str(vstretch)])
    end

    xlabel([num2str(param) '  ' int2str(hc) '  ' int2str(n)])

    figure
    plot(s,dz,'+')
    xlabel('s')
    ylabel('\Deltaz (m)')
end

