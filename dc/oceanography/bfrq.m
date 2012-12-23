function [N2] = bfrq(s,t,z,lat)
    % Calculate buoyancy frequency @ midpts of input grid z
    % Input
    %           s -> salinity (psu)
    %           t -> temperature (celsius)
    %           z -> depths
    %           lat -> latitude

    lz = length(z);
	N2 = (-9.81./sw_smow(t(1:lz-1))).*ones(lz-1,1);
	pref = sw_pres((z(1:lz-1)+z(2:lz))/2,lat);
    p = sw_pres(z,lat);
    
	for k=1:lz-1
        N2(k) = N2(k)*(sw_pden(s(k),t(k),p(k),pref(k))-sw_pden(s(k+1),t(k+1),p(k+1),pref(k)))/(z(k+1)-z(k));
	end
    
    % Account for missing values in T at large depths
    ind = sort(find(isnan(N2) == 1));
    if isempty(ind)
        return
    else
        fprintf('\n NaN detected in N2.');
        %if ind(1) >= lz-2
            fprintf(' Filling %d values with nearest value @ depth %d', length(ind),ind(1)-1);
            ind
            N2(ind) = N2(ind(1)-1);
            % end
    end
    %end