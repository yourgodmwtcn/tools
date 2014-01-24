function xg = resolution2grid(x,resol,r);

% xg = resolution2grid(x,resol,r)
%
% works out grid coordinates along one axis of a cartesian grid that allow the
% grid to have specified resolution at specified points resol(x), with a max
% ratio between adjacent cell widths of r.
%
% neil banas, may 2011

do_plot = 0;

x_orig = x;
resol_orig = resol;

dxg = [];
% step through segments, assembling the grid as we go
for i=1:length(x)-1
	if resol(i) == resol(i+1)
		N = ceil((x(i+1) - x(i)) / resol(i)); % number of points that span the interval
		D = N * resol(i); % what the interval distance has to be for N to be an integer
		x(i+1) = x(i) + D; % adjust the far end of the interval to make this work
		dxg = [dxg; resol(i) .* ones(N,1)]; % add N cells to dxg 
	else
		narrowing = (resol(i) > resol(i+1)); % as opposed to widening (boolean)
		dmax = max(resol(i), resol(i+1));
		dmin = min(resol(i), resol(i+1));
		NT = ceil((log(dmax) - log(dmin))/log(r) - 1); % number of points in the transition zone
		reff = exp((log(dmax) - log(dmin)) / (NT+1)); % r adjusted to make NT an integer
		DT = (dmax - reff * dmin) / (reff-1); % distance required for the transition zone
		DF = max(x(i+1) - x(i) - DT, 0); % remaining distance: the flat part at dmax
		if (narrowing)
			NF = floor(DF/dmax); % number of points in the flat part
		else
			NF = ceil(DF/dmax);
		end
		if i == length(x) - 1 & ~narrowing
			dmax = DF/NF; % in the last interval, if widening, adjust resolution rather than spacing
			resol(i+1) = dmax;
		else
			DF = NF * dmax; % make sure DF is an integer number of cells
			x(i+1) = x(i) + DF + DT; % adjust the far end of the interval to make everything fit
		end
		if (narrowing)
			dxg = [dxg; dmax .* ones(NF,1)];
			for j=1:NT
				dxg = [dxg; dmax/reff^j];
			end
		else
			for j=1:NT
				dxg = [dxg; dmin*reff^j];
			end
			dxg = [dxg; dmax .* ones(NF,1)];
		end
	end
end
xg = x(1) + [0; cumsum(dxg)];

if do_plot
	figure
	plot(x,resol,'ro');
	hold on
	for i=1:length(xg)
		plot([xg(i) xg(i)],ylim,'color',[0.8 0.8 0.8]);
	end
	plot(0.5.*(xg(1:end-1)+xg(2:end)),dxg,'k.');
	plot(x_orig,resol_orig,'r*');
	plot(x,resol,'ro');
end

%disp(length(xg))