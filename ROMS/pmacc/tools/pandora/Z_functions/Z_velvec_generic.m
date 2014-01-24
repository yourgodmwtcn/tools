function [] = Z_velvec_generic(lon,lat,u,v,poschar,utyp)
% 11/30/2011 Parker MacCready
%
% This adds surface velocity vectors to an existing plot.  It takes care to
% account for scaling of Cartesian lat/lon grids, so that a current to the
% NE actually points with an angle of 45 degrees on the plot.
%
% This requires input of the coordinates and velocity, and is best applied
% to regularly-spaced data
%
% poschar is the location of the SCALE ARROWS
%   'ul','ur','ll','lr' for the corners, or
%   'none' for none
%
% utyp = typical velocity

hold on

aa = axis; Dlat = aa(4)-aa(3); Dlon = aa(2)-aa(1);
darscale = 1/cos(pi*mean(aa(3:4))/180);
if nargin==5; utyp = 1; end;
ufact = Dlat/(30*utyp); % automatic scaling of arrow size;
quiver(lon,lat,ufact*u*darscale,ufact*v,0,'k');

% add veocity scale arrows
uscale = utyp; % the velocity of the scale arrows
switch poschar
    case 'none'
        % do nothing
    otherwise
        [xt,yt] = Z_lab(poschar);
end
switch poschar
    case {'ul'}
        deltax = 0; deltay = -Dlat/20;
        quiver(xt+deltax,yt,ufact*uscale*darscale,ufact*0,0,'k');
        quiver(xt+deltax,yt,ufact*0*darscale,ufact*uscale,0,'k');
        text(xt+deltax,yt+deltay,[num2str(uscale),' m s^{-1}'],'color','k');
    case {'ur'}
        deltax = -Dlon/4; deltay = -Dlat/20;
        quiver(xt+deltax,yt,ufact*uscale*darscale,ufact*0,0,'k');
        quiver(xt+deltax,yt,ufact*0*darscale,ufact*uscale,0,'k');
        text(xt+deltax,yt+deltay,[num2str(uscale),' m s^{-1}'],'color','k');
    case {'ll'}
        deltax = 0; deltay = Dlat/20;
        quiver(xt+deltax,yt+deltay,ufact*uscale*darscale,ufact*0,0,'k');
        quiver(xt+deltax,yt+deltay,ufact*0*darscale,ufact*uscale,0,'k');
        text(xt+deltax,yt,[num2str(uscale),' m s^{-1}'],'color','k');
    case {'lr'}
        deltax = -Dlon/4; deltay = Dlat/20;
        quiver(xt+deltax,yt+deltay,ufact*uscale*darscale,ufact*0,0,'k');
        quiver(xt+deltax,yt+deltay,ufact*0*darscale,ufact*uscale,0,'k');
        text(xt+deltax,yt,[num2str(uscale),' m s^{-1}'],'color','k');
    case 'none'
        % do nothing
end
