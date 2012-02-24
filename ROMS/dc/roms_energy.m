% function [EKE,MKE,PE] = roms_energy(fname,tindices)

% CHANGELOG
% First actual version                                          23 Feb 2012

% Todo list
% 1) add netcdf output
% 2) Is there a bug in PE calculation?
% 3) Try lower values of slab to see which works best
% 4) Add w contribution to KE

function [EKE,MKE,PE] = roms_energy(fname,tindices)

% input
fname = 'his';
fname = find_file(fname);
tindices = [1 Inf];

% parameters
vinfo = ncinfo(fname,'u');
dim   = length(vinfo.Size); 
slab  = 40;

warning off
grid = roms_get_grid(fname,fname,0,1);
warning on

%% parse input

% FROM mod_movie.m - propagate changes back
if ~exist('tindices','var') || isempty(tindices)
    tindices = [1 Inf];
    dt = 1;
else
    switch length(tindices)
        case 1
            dt = 1;
            tindices(2) = tindices(1);

        case 2
            dt = 1;

        case 3
            dt = tindices(2);
            tindices(2) = tindices(3);
            tindices(3) = NaN;
    end
end

if tindices(2) < tindices(1)
    tindices(2) = tindices(1)+tindices(2);
end
if isinf(tindices(2)), tindices(2) = vinfo.Size(end); end

stride = [1 1 1 dt];

if (tindices(2)-tindices(1)) == 0
    iend = 1;
    dt = tindices(2);
else
    iend   = ceil((tindices(2)-tindices(1))/slab/dt);
end
% END mod_movie section
nt = ceil(tindices(2)-tindices(1)+1)/dt;

%% read data

% caps indicates domain integrated values
EKE = nan(nt,1);
MKE = EKE;
PE  = EKE;

R0  = ncread(fname,'R0');
time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;

try
    cpb = progressbar();
catch ME
    cpb = [];
end

for i=0:iend-1
    % FROM mod_movie.m - propagate changes back
    % start and count arrays for ncread : corrected to account for stride
    read_start = ones(1,dim);
    read_count = Inf(1,dim);

    if i == (iend-1)
        read_count(end) = ceil((tindices(2)-slab*(i))/dt);
    else
        read_count(end) = ceil(slab/dt);%ceil(slab*(i+1)/dt);
    end

    if i == 0
        read_start(end) = tindices(1);
    else
        read_start(end) = slab*i + 1;
    end

    if (iend-1) == 0, read_count(end) = ceil((tindices(2)-tindices(1))/dt)+1; end 
    % END mod_movie section
    
    if isempty(cpb), fprintf('\nReading Data...\n'); end
    u   = ncread(fname,'u',read_start,read_count,stride); pbar(cpb,i+1,1,iend,4);
    v   = ncread(fname,'v',read_start,read_count,stride); pbar(cpb,i+1,2,iend,4);
    w   = ncread(fname,'w',read_start,read_count,stride); pbar(cpb,i+1,3,iend,4);
    rho = R0 + ncread(fname,'rho',read_start,read_count,stride); pbar(cpb,i+1,4,iend,4);
	if isempty(cpb), fprintf('\n Done reading data... \n'); end
    
    % mean fields
    um = mean(u,2);
    vm = mean(v,2);
    wm = mean(w,2);
    rm = mean(rho,2);

    % eddy fields
    up = sub_along_mean(u);
    vp = sub_along_mean(v);
    wp = sub_along_mean(w);
    rp = sub_along_mean(rho);
    
    eke = 0.5*rho(1:end-1,1:end-1,:,:).*(up(:,1:end-1,:,:).^2 + vp(1:end-1,:,:,:).^2);
    mke = 0.5*rm(1:end-1,:,:,:).*(um(:,:,:,:).^2 + vm(1:end-1,:,:,:).^2);
    pe = 9.81*rho.*repmat(permute(grid.z_r,[3 2 1]),[1 1 1 size(rho,4)]);
    
    s = size(u);
    
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end);

    EKE(tstart:tend-1) = squeeze(sum(sum(sum(eke,1),2),3));%zeros(size(eke,4),1);
    MKE(tstart:tend-1) = squeeze(sum(sum(sum(mke,1),2),3));%zeros(size(mke,4),1);
    PE(tstart:tend-1)  = squeeze(sum(sum(sum(pe,1),2),3));%zeros(size(pe,4),1);
    
end

if ~isempty(cpb)
    cpb.stop();
end

%% plot
figure;
plot(time,PE/nanmax(PE(:)),'b'); hold on;
plot(time,EKE/nanmax(EKE(:)),'r');
plot(time,MKE/nanmax(MKE(:)),'k');
ylabel('Energy');
xlabel('Time (days)');
legend('PE','EKE','MKE');

figure;
plot(time,PE);
ylabel('Energy');
xlabel('Time (days)');
legend('PE');

%% old code
%     % vectorize this?
%     for i=1:s(4)
%         temp = eke(:,:,:,i);
%         int_eke(i) = sum(temp(:));
%         temp = pe(:,:,:,i);
%         int_pe(i) = sum(temp(:));
%     end  
%     

%     % eddy terms
%     up = u - repmat(mean_u,[1 1 1 s(4)]);
%     vp = v - repmat(mean_v,[1 1 1 s(4)]);

%% local functions
function [ap] = sub_along_mean(a)
    ap = bsxfun(@minus,a,mean(a,2));
    
function [] = pbar(cpb,i,j,imax,jmax)
    if ~isempty(cpb)
        txt = sprintf(' Progress: i=%d, j=%d',i,j);
        progressbarupdate(cpb,(jmax*(i-1)+j)/(imax*jmax)*100,txt);
    end