function [eke,tke,tpe,pep,zpe,erf,mrf,ntt,rbbb,tavg] = engflcal(fname,tlow,thigh,xlow,xhigh,zzero,rbbb,nts)

%   compute energy for time window (tlow,thigh), x window (xlow,xhigh)
%       also does heat flux at centrer
%   call as
%       [eke,tke,tpe,epe,zpe,ntt,rbbb,tavg] = engflcal(fname,tlow,thigh,xlow,xhigh,zzero,rbbb,nts)

%   6/16/2011      8/17/2011
%   corrected for v grid 8/29/2011

g = 9.8;
rz = 1027;
tzero = 24.0;
bet = 1.7e-4;
eps = 0.0001;



ncid = netcdf.open(fname,'NC_NOWRITE');
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'))/86400;

if time(1) > eps
    iadd = 0;
else
    iadd = 1;
end
%time = squeeze(time((1+iadd):nts:end));

ntime = length(time);
%[tt, itl] = min(abs(time-tlow));
%[tt,ith] = min(abs(time-thigh));
%itt = itl:ith;
itt = find(time > (tlow+eps) & time < (thigh+eps));
itl = min(itt);
ith = max(itt);
ntt = length(itt);
tavg = mean(time(itt));


time = squeeze(time((1+iadd):nts:end));


xr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_rho'))/1000.;
[nxr,nyr] = size(xr);
xr = xr(:,5);
%[xx,irho] = min(abs(xr-xc));
dx =(xr(2) - xr(1));
ixrl = round(1.5 + xlow/dx);
ixrh = round(1.5 + xhigh/dx);
nxrr = length(ixrl:ixrh);
icent = round(mean(1:nxrr));

dx = dx*1000;


srho = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'s_rho'));
nz = length(srho);
csr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cs_r'));
hc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'));
theta_s = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'));
theta_b = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'));
Tcline = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Tcline'));
vxform = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform'));
vstretch = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));

h = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'),[0,5],[nxr,1]);
h = squeeze(h);
hh = h(ixrl:ixrh);

zz = NaN*ones(nxrr,nz);
dz = zz;

for ix = 1:nxrr
    %zz(ix,:) =  hc*srho + (hh(ix)-hc)*csr;
    if round(vxform) ==1
        zz(ix,:) =  hc*srho + (hh(ix)-hc)*csr;
    else
        sss = (hc*srho + hh(ix)*csr)/(hc+hh(ix));
        zz(ix,:) = hh(ix)*sss;
    end
    dz(ix,2:end-1) = (zz(ix,3:end) - zz(ix,1:end-2))/2;
    dz(ix,1) = (zz(ix,2) + zz(ix,1))/2 + hh(ix);
    dz(ix,end) = -(zz(ix,end) + zz(ix,end-1))/2;
end




vi = NaN*ones(nxrr,nyr-2,nz,ntt);
%itl-1+iadd
%ntt
v = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[ixrl-1,0,0,itl-1],[nxrr,nyr-1,nz,ntt]);
%size(v)
v = squeeze(v);
%size(v)
%v = (v(:,:,:,(1+iadd):nts:end));
v(:,nyr-1,:,:) = v(:,1,:,:);
vi(:,1:nyr-2,:,:) = 0.5*(v(:,1:nyr-2,:,:) + v(:,2:nyr-1,:,:));
clear v
vi = squeeze(vi(:,:,:,1:nts:end));

[njj1,njj2,njj3,ntt] = size(vi);

vimt = mean(vi,4);          % time mean
vimt = squeeze(vimt);
vimty = mean(vimt,2);       % t,y mean
clear vimt
vimty = squeeze(vimty);
ery = 0*ones(nyr-2,nz);
ermy = ery;

%vp = vi*NaN;
vsq = 0*vimty;
vtke = vsq;


for ix = 1:nxrr
    for iz = 1:nz
        vm = vimty(ix,iz);
        da = dz(ix,iz)*dx;
        for iy = 1:nyr-2
            for it = 1:ntt  
                vsq(ix,iz) = vsq(ix,iz) + da*(vi(ix,iy,iz,it)-vm)^2;
                vtke(ix,iz) = vtke(ix,iz) + da*(vi(ix,iy,iz,it))^2;
            end
        end
    end
end


nnv = ntt*(nyr-2);
vsq = vsq/nnv;
vtke = vtke/nnv;

clear vi
clear vimty




ui = NaN*ones(nxrr,nyr-2,nz,ntt);

u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[ixrl-2,1,0,itl-1],[nxrr+1,nyr-2,nz,ntt]);
u = squeeze(u);
%u = (u(:,:,:,(1+iadd):nts:end));
ui(1:nxrr,:,:,:) = 0.5*(u(1:nxrr,:,:,:) + u(2:nxrr+1,:,:,:));
clear u
ui = squeeze(ui);
ui = squeeze(ui(:,:,:,1:nts:end));




%if ntt+ 0.001 > 1 
uimt = mean(ui,4);
uimt = squeeze(uimt);
uimty = mean(uimt,2);
clear uimt
uimty = squeeze(uimty);



usq = 0*uimty;
utke = usq;

for ix = 1:nxrr
    for iz = 1:nz
        um = uimty(ix,iz);
        da = dz(ix,iz)*dx;
        for iy = 1:nyr-2
            for it = 1:ntt
                
                usq(ix,iz) = usq(ix,iz) + da*(ui(ix,iy,iz,it)-um)^2;
                utke(ix,iz) = utke(ix,iz) + da*(ui(ix,iy,iz,it))^2;
            end
        end
    end
end

ucent = squeeze(ui(icent,:,:,:));

clear ui
clear uimty

nnu = (nyr-2)*ntt;
usq = usq/nnu;
utke = utke/nnu;


a = 0;
for ix = 1:nxrr
    for iz = 1:nz;
        a = a + dz(ix,iz)*dx;
    end
end


%%%
%  Now do PE calculations

ifirst = 0;
if 0*rbbb ~= 0
    rbbb = 0;
    ifirst = 1;
end


temp = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),[ixrl-1,1,0,itl-1],[nxrr,nyr-2,nz,ntt]);
temp = squeeze(temp);
rho = -bet*rz*(temp-tzero) -rbbb;
rho = squeeze(rho(:,:,:,1:nts:end));

%if (ntt+0.0001) > 1 
 rmt = mean(rho,4);
 rmt = squeeze(rmt);
 rmty = mean(rmt,2);
 rmty = squeeze(rmty);
 


if ifirst == 1
   
    for ix = 1:nxrr
        for iz = 1:nz
            rbbb = rbbb + rmty(ix,iz)*dz(ix,iz)*dx;
        end
    end
    
    rbbb = rbbb/a;
    rho = -bet*rz*(temp-tzero) -rbbb;
    rmt = mean(rho,4);
    rmt = squeeze(rmt);
    rmty = mean(rmt,2);
    rmty = squeeze(rmty);
    disp('Time(days) eke tke,pep,tpe, zpe')
end

clear temp


tpex = 0*rmty;
%mpex = tpex;
pepx = tpex;

for ix = 1:nxrr
    for iz = 1:nz
        rm = rmty(ix,iz);
        da = dz(ix,iz)*dx;
        %tpex(ix,iz) = da*(zz(ix,iz)-zzero)*rm;
        for iy = 1:nyr-2
            for it = 1:ntt
                pepx(ix,iz) = pepx(ix,iz) + (da*(zz(ix,iz)-zzero)*rho(ix,iy,iz,it))^2;
                tpex(ix,iz) = tpex(ix,iz) + da*(zz(ix,iz)-zzero)*rho(ix,iy,iz,it);
                %mpex(ix,iz) = mpex(ix,iz) + da*(zz(ix,iz)-zzero)*(rm);
            end
        end
    end
end

%mpex = mpex/(nnu);
tpex = tpex/(nnu);
pepx = real(sqrt(pepx/nnu  -tpex.*tpex));
tpe = g*sum(sum(tpex))/(a*rz);
%mpe = g*sum(sum(mpex))/(a*rz);
pep = g*sum(sum(pepx))/(a*rz);

rhocent = squeeze(rho(icent,:,:,:));

%   compute density flux

clear rho
clear rmt
clear rmty

ucentm = mean(ucent,3);         % time mean
rcentm = mean(rhocent,3);
ucentmty = squeeze(mean(ucentm,1));      % y,t mean
rcentmty = squeeze(mean(rcentm,1));
hcent = sum(dz(icent,:));       % water depth


for it = 1:ntt
    for iz = 1:nz
        ery(:,iz) = ery(:,iz) + dz(icent,iz)*ucent(:,iz,it).*rhocent(:,iz,it);
        ermy(:,iz) =ermy(:,iz) + dz(icent,iz)*ucentmty(iz).*rcentmty(iz);
    end
end
ery = squeeze(mean(ery,1));                  % y mean
ermy = squeeze(mean(ermy,1));

erf = sum(ery)/(hcent*ntt);        % total flux
mrf = sum(ermy)/(ntt*hcent);       % flux due to y,t mean

clear rhocent
clear ucent


%%%%%%%%%%%%%
% Now do zeta PE


zet = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zeta'),[ixrl-1,1,itl-1],[nxrr,nyr-2,ntt]);
zet = squeeze(zet);
zet = squeeze(zet(:,:,1:nts:end));

zetmt = mean(zet,3);
zetmt = squeeze(zetmt);
zetmty = mean(zetmt,2);
zetmty = squeeze(zetmty);

zetpex = 0*ones(nxrr,1);
%zetpew = zetpex;

for ix = 1:nxrr
   
        
        for iy = 1:nyr-2
            for it = 1:ntt
                %zetpex(ix,iz) = zetpex(ix,iz) + (da*(zz(ix,iz)-zzero)*rho(ix,iy,iz,it))^2;
                zetpex(ix) = zetpex(ix) + zet(ix,iy,it)^2;
                %zetpew(ix) = zetpew(ix) -zzero*zet(ix,iy,it);
                %mpex(ix,iz) = mpex(ix,iz) + da*(zz(ix,iz)-zzero)*(rm);
            end
        end
    
end

clear zet

zetpex = zetpex/(ntt*(nyr-2));
zpe = sum(zetpex)/nxrr;
%zpew = sum(zetpew)/nxrr;
xspan = dx*nxrr;
hbar = a/xspan;
zpe = zpe*g/(2*hbar) ;


netcdf.close(ncid)
tmm = (tlow + thigh)/2;

%nn = ntt*(nyr-3);
%ekez = 0.5*(vsq+usq)/(nn);     
eke = 0.5*(sum(sum(vsq + usq)))/(a);
tke = 0.5*(sum(sum(vtke + utke)))/(a);

[tavg eke, tke pep tpe zpe]

