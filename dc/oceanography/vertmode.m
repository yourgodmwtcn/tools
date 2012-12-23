function [Vmode, Hmode, c] = vertmode(N2, Z, n, make_plot)
	
	% function [Vmode, Hmode] = vertmode(N2, Z, n)
	% Takes input N2 -> Buoyancy frequency squared (@ mid pts of Z)
	%              Z -> Vertical co-ordinate
	%              n -> No. of modes to isolate
    %                   Max. number of nodes is length(Z)-1
	% Returns
	%             Vmode -> Vertical structure for vertical velocity
	%             Hmode -> Vertical structure for horizontal velocities, pressure
    %              c(i) -> Gravity Wave speed of i-th mode 
    
    if ~exist('make_plot','var'), make_plot = 1; end
    
	lz = length(Z);
    
    if n > (lz-1)
        fprintf('\n n too big. Reducing to (length(Z) - 1)');
        n = lz - 1;
    end

	% Following Chelton (1998)
	% Q'AW = (lambda)W
	
	D(1) = Z(2)/2;
	D(2:lz-1) = (Z(3:lz)-Z(1:lz-2))/2;
	D(lz) = (Z(lz) - Z(lz-1))/2;
	Zmid = (Z(1:lz-1) + Z(2:lz))/2;
    
    Q1 = zeros(lz-1,1);
	
    for k = 1:lz-1
        Q1(k) = 2/(D(k)*D(k+1)*(D(k) + D(k+1)));
    end
	
	Q = diag(Q1./N2);
	A = -1*(diag(D(1:lz-1)) + diag(D(2:lz))) + diag(D(1:lz-2),1) + diag(D(3:lz),-1);
	
	[G,e] = eig(Q*A);
	
	c = sqrt(-1./diag(e));
	[~,ind] = sort(c); % ascending order
    
    F=zeros(lz-2,lz-1);
	for i=1:lz-1
        F(:,i) = (c(i)^2/9.81)*(diff(G(:,i))./diff(Zmid));
	end
    
    Hmode = fliplr(F(:,ind(lz-n:lz-1)));
    Vmode = fliplr(G(:,ind(lz-n:lz-1)));
    
    % Fill in Hmode at lowest depth
    Hmode(lz-1,:) = Hmode(lz-2,:);
    c = flipud(c(ind(lz-n:lz-1)));
    
    % Normalize by max. amplitude
    %Hmode = Hmode./repmat(max(abs(Hmode)),length(Hmode),1);
    %Vmode = Vmode./repmat(max(abs(Vmode)),length(Vmode),1);
    
    % Normalize by energy (following Wunsch(1999)
    norm = sqrt(sum(avg1(Hmode).^2.*repmat(diff(Zmid),1,n)));
    Vnorm = Vmode./repmat(norm,lz-1,1);
    Hnorm = Hmode./repmat(norm,lz-1,1);

    % check normalization
    hchk = sum(avg1(Hnorm).^2.*repmat(diff(Zmid),1,n));
    vchk = sum(avg1(Vnorm.*repmat(N2,1,n)).^2.*repmat(diff(Zmid),1,n));
    
    fprintf('\n u mode normalization: \n');
    disp(hchk);
    fprintf('\n w mode normalization: \n');
    disp(vchk);
    
    Hmode = Hnorm;
    Vmode = Vnorm;
    
    % Plot first n modes
    if make_plot
        figure;
        subplot(121)
        plot(Vmode,Zmid);
        set(gca,'ydir','reverse');
        hold on;
        %xlim([-1.5 1.5]);
        %ylim([0 600]);
        plot([0 0],ylim, 'k-');
        plot(xlim,[4000 4000], 'k--');
        ylabel('Z (m)');
        title('Vertical Structure of Vertical Velocity');
        legend(num2str((1:n)'));
        beautify;

        subplot(122)
        plot(Hmode,Zmid);
        %xlim([-1.5 1.5]);
        %ylim([0 600]);
        set(gca,'ydir','reverse');
        hold on;
        plot([0 0],ylim, 'k-');
        plot(xlim,[4000 4000], 'k--');
        ylabel('Z (m)');
        title('Vertical Structure of Horizontal Velocity');
        legend(num2str([1:n]'));
        beautify;
    end