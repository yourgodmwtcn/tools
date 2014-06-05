% compute baroclinicity - better defined as metric of vertical
% uniformity of a profile.
%     bc = (KE_total - KE_depthavg)/KE_tot 
%          where KE = int_z 0.5 * profile.^2 dz

function [bc] = baroclinicity(zvec, profile)

    pmean = nanmean(profile);

    ketot = abs(0.5 * trapz(zvec, profile.^2));
    keda = abs(0.5 * trapz(zvec, pmean^2 .* ones(size(profile))));

    bc = (ketot - keda)/ketot;