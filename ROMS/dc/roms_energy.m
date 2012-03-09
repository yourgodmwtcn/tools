% function [EKE,MKE,PE] = roms_energy(fname,tindices,ntavg,mean_index)
% Use mean_index to say which dirn. you want to take the mean for defining
% mean, eddy contributions.
% Normalizes energy by horizontal area

% CHANGELOG
% First actual version                                          23 Feb 2012

% Todo list
% 1) add netcdf output
% 2) Is there a bug in PE calculation?
% 3) Try lower values of slab to see which works best
% 4) Add w contribution to KE

function [EKE,MKE,PE] = roms_energy(fname,tindices,ntavg,mean_index)

if ~exist('tindices','var'), tindices = [1 Inf]; end
if ~exist('mean_index','var'), mean_index = 2; end
if ~exist('ntavg','var'), ntavg = 3; end % average over ntavg timesteps

% parameters
vinfo = ncinfo(fname,'u');
dim   = length(vinfo.Size); 
slab  = 60;

warning off
grid = roms_get_grid(fname,fname,0,1);
warning on

area = max(grid.x_rho(:)-grid.x_rho(1))*max(grid.y_rho(:)-grid.y_rho(1));

% parse input
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

%% read data

% caps indicates domain integrated values
EKE = nan(floor(nt/ntavg),1);
MKE = EKE;
PE  = EKE;

R0  = ncread(fname,'R0');
time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;

zrho = permute(grid.z_r,[3 2 1]);
try
    cpb = progressbar();
catch ME
    cpb = [];
end

% create output file
% outname = 'ocean_der.nc';
% xname = 'x_pv'; yname = 'y_pv'; zname = 'z_pv'; tname = 'ocean_time';
% try
%     nccreate(outname,'pv','Dimensions', {xname s(1)-1 yname s(2)-2 zname s(3)-1 tname length(tpv)});
%     nccreate(outname,xname,'Dimensions',{xname s(1)-1});
%     nccreate(outname,yname,'Dimensions',{yname s(2)-2});
%     nccreate(outname,zname,'Dimensions',{zname s(3)-1});
%     nccreate(outname,tname,'Dimensions',{tname length(tpv)});
%     
%     ncwriteatt(outname,'pv','Description','Ertel PV calculated from ROMS output');
%     ncwriteatt(outname,'pv','coordinates','x_pv y_pv z_pv ocean_time');
%     ncwriteatt(outname,'pv','units','N/A');
%     ncwriteatt(outname,xname,'units',ncreadatt(fname,'x_u','units'));
%     ncwriteatt(outname,yname,'units',ncreadatt(fname,'y_u','units'));
%     ncwriteatt(outname,zname,'units','m');
%     ncwriteatt(outname,tname,'units','s');
%     fprintf('\n Created file : %s\n', outname);
% catch ME
%     fprintf('\n Appending to existing file.\n');
% 

for i=0:iend-1
    % FROM mod_movie.m - propagate changes back
    % start and count arrays for ncread : corrected to account for stride
    
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    
    if isempty(cpb), fprintf('\nReading Data...\n'); end
    u   = ncread(fname,'u',read_start,read_count,stride); pbar(cpb,i+1,1,iend,4);
    v   = ncread(fname,'v',read_start,read_count,stride); pbar(cpb,i+1,2,iend,4);
    %w   = ncread(fname,'w',read_start,read_count,stride); pbar(cpb,i+1,3,iend,4);
    rho = R0 + ncread(fname,'rho',read_start,read_count,stride); pbar(cpb,i+1,4,iend,4);
	if isempty(cpb), fprintf('\n Done reading data... \n'); end
    
    % mean fields - average over 4 timesteps
    um = time_mean(u,ntavg,mean_index);
    vm = time_mean(v,ntavg,mean_index);
    %wm = mean(w,2);
    rm = time_mean(rho,ntavg,mean_index);
    
    s = size(u);

    % eddy fields
    ind1 = ceil(ntavg/2)  :ntavg:s(4)-mod(s(4),ntavg);
    
    if ntavg == 1
        ind2 = ind1; 
    else
        ind2 = ceil(ntavg/2+1) :ntavg:s(4)-mod(s(4),ntavg);
    end
    
    % pull out rho at timesteps where i'm calculating eddy fields.
    rho  = (rho(:,:,:,ind1) + rho(:,:,:,ind2))/2;
    
    up = bsxfun(@minus,(u(:,:,:,ind1) + u(:,:,:,ind2))/2,um);
    vp = bsxfun(@minus,(v(:,:,:,ind1) + v(:,:,:,ind2))/2,vm);
    %wp = bsxfun(@minus,w,wm);
    rp = bsxfun(@minus,rho,rm);
    
    % average so that everything lands up on interior-rho points
    up = (up(1:end-1,2:end-1,:,:) + up(2:end,2:end-1,:,:))/2;    
    vp = (vp(2:end-1,1:end-1,:,:) + vp(2:end-1,2:end,:,:))/2;
    if mean_index == 1
        um = um(:,2:end-1,:,:);
        vm = (vm(:,1:end-1,:,:) + vm(:,2:end,:,:))/2;
    elseif mean_index == 2
        um = (um(1:end-1,:,:,:) + um(2:end,:,:,:))/2;
        vm = vm(2:end-1,:,:,:);
    end
    
    % now calculate energy terms
    eke = 0.5*rho(2:end-1,2:end-1,:,:).*(up.^2 + vp.^2)./area; % SLOW?!
    mke = 0.5*bsxfun(@times,rho(2:end-1,2:end-1,:,:),(um.^2 + vm.^2))./area;
    %oke = rho(2:end-1,2:end-1,:,:).*(bsxfun(@times,up,um)+ bsxfun(@times,vp,vm))./area;
    pe  = 9.81*bsxfun(@times,rho,zrho)./area;
    
    tstart = ceil(read_start(end)/ntavg);
    tend   = floor(tstart + s(4)/ntavg -1);

    time_E(tstart:tend,1) = (time(read_start(end)+ind1-1) + time(read_start(end)+ind2-1))/2;
    EKE(tstart:tend) = domain_integrate(eke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
    MKE(tstart:tend) = domain_integrate(mke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
   % OKE(tstart:tend) = domain_integrate(oke,grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));
    PE(tstart:tend)  = domain_integrate(pe(2:end-1,2:end-1,:,:),grid.x_rho(1,2:end-1)',grid.y_rho(2:end-1,1),grid.z_r(:,1,1));  
    
    % Write to netcdf file here
end

if ~isempty(cpb)
    cpb.stop();
end

%% Calculate growth rate
k=1;
jump = 5; % fit 5 consecutive points

for i=1:1:length(EKE)-jump
    A(k,:) = polyfit(time_E(i:i+jump),log(EKE(i:i+jump)),1);
    time_A(k,:) = (time_E(i+jump)+time_E(i))/2;
    k=k+1;
end

figure
subplot(211)
plot(time_A,A(:,1),'b*-')
liney(0)
ylabel('Growth Rate (d^{-1})')
xlabel('Time (days)');

% Verify
eke2 = exp(A(:,1).*time_A + A(:,2));
subplot(212)
plot(time_E,(EKE),'b*');
hold on
plot(time_A,eke2,'r');
ylabel('Energy');
xlabel('Time (days)');
title('Verification');
legend('Original','Fit');

A = A(:,1);

%% plot

figure;
hold on;
plot(time_E,EKE,'r');
plot(time_E,MKE,'k');
%plot(time,OKE,'m');
ylabel('Energy');
xlabel('Time (days)');
legend('EKE','MKE');

figure;
plot(time_E,PE);
ylabel('Energy');
xlabel('Time (days)');
legend('PE');

% write to file
ax = 'xyzt';
fname = ['energy-avg-' ax(mean_index) '.mat']; 
save(fname,'time_E','PE','EKE','MKE','A','time_A','ntavg');

%% local functions

function [datam] = time_mean(data,n,mean_index)
    for ii = 1:n:size(data,4)-n+1
        ind = ceil(ii/n);
        datam(:,:,:,ind) = mean(mean(data(:,:,:,ii:ii+n-1),4),mean_index);
    end
        
function [] = pbar(cpb,i,j,imax,jmax)
    if ~isempty(cpb)
        txt = sprintf(' Progress: i=%d, j=%d',i,j);
        progressbarupdate(cpb,(jmax*(i-1)+j)/(imax*jmax)*100,txt);
    end