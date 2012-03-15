function [encontot,enconm,cb,cs] = enconcal(fname,tlow,thigh,xlow,xhigh,rbbb,nts);

%   compute the energy conversions gfor the given time, space window
%   call from engtsflc.m as
%       [encontot,enconm,cb,cs] = enconcal(fname,tlow,thigh.xlow,xhigh,rbbb,nts);
%   where
%       fname = file name (string)
%       tlow, thigh = time window
%       xlow , xhigh = x window
%       rbbb = reference density
%       nts = 1  stride

%       9/22/2011

g = 9.8;
rz = 1027;
tzero = 14.0;
bet = 1.7e-4;
eps = 0.0001;

const = -g/rz;

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

xr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_rho'))/1000.;
[nxr,nyr] = size(xr);
xr = xr(:,5);
%[xx,irho] = min(abs(xr-xc));
dx =(xr(2) - xr(1));
ixrl = round(1.5 + xlow/dx);
ixrh = round(1.5 + xhigh/dx);
nxrr = length(ixrl:ixrh);
dx = dx*1000;



srho = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'s_rho'));
ds = srho(2) - srho(1);
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
hx = 0*hh;
hx(2:end-1) = (hh(3:end) - hh(1:end - 2))/dx;


zz = NaN*ones(nxrr,nz);
dz = zz;
zs = NaN*ones(nz,1);

for ix = 1:nxrr
    %zz(ix,:) =  hc*srho + (hh(ix)-hc)*csr;
    if round(vxform) ==1
        zz(ix,:) =  hc*srho + (hh(ix)-hc)*csr;
    else
        sss = (hc*srho + hh(ix)*csr)/(hc+hh(ix));
        zz(ix,:) = hh(ix)*sss;
    end
     zs(2:end-1) = 0.5*(zz(ix,3:end) - zz(ix,1:end-2))/ds; 
    zs(1) = (0.5*(zz(ix,2)+zz(ix,1)) +hh(ix))/ds;
    zs(end) = -0.5*(zz(ix,end-1)+zz(ix,end))/ds;
    dxc(ix,:) = -hx(ix)*csr./zs;
    
    dz(ix,2:end-1) = (zz(ix,3:end) - zz(ix,1:end-2))/2;
    dz(ix,1) = (zz(ix,2) + zz(ix,1))/2 + hh(ix);
    dz(ix,end) = -(zz(ix,end) + zz(ix,end-1))/2;
end


temp = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),[ixrl-1,1,0,itl-1],[nxrr,nyr-2,nz,ntt]);
temp = squeeze(temp);
rho = -bet*rz*(temp-tzero) -rbbb;
clear temp
rho = squeeze(rho(:,:,:,1:nts:end));

wread = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'w'),[ixrl-1,1,0,itl-1],[nxrr,nyr-2,nz+1,ntt]);
wread = squeeze(wread);
wi = 0.5*(wread(:,:,1:end-1,:) + wread(:,:,2:end,:));
clear wread
wi = squeeze(wi(:,:,:,1:nts:end));


wmean = 0*ones(nxrr,nz);
rmean = wmean;

if ntt*0.98 > 1
    wmeant = mean(wi,4);
    wmean = mean(wmeant,2); % y,t average
    rmeant = mean(rho,4);
    rmean = mean(rmeant,2);
else
    wmean = squeeze(mean(wi,2));     
    rmean = squeeze(mean(rho,2));
end

cmean = 0*ones(1,nz);
ctot = cmean;
area = 0;


for iz = 1:nz
    area = area + mean(dz(:,iz));
    cmean(iz) = mean(dz(:,iz).*rmean(:,iz).*wmean(:,iz));
    if ntt*0.98 > 1
        flt = mean(wi(:,:,iz,:).*rho(:,:,iz,:),4);
                
    else
        flt = wi(:,:,iz).*rho(:,:,iz);
    end
    cmeanty = mean(flt,2);
    ctot(iz) = mean(dz(:,iz).*cmeanty);
end


enconm = const*sum(cmean)/area;
encontot = const*sum(ctot)/area;
clear rho



ui = NaN*ones(nxrr,nyr-2,nz,ntt);
u = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[ixrl-2,1,0,itl-1],[nxrr+1,nyr-2,nz,ntt]);
u = squeeze(u);
ui(1:nxrr,:,:,:) = 0.5*(u(1:nxrr,:,:,:) + u(2:nxrr+1,:,:,:));
clear u
ui = squeeze(ui);

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
vi = squeeze(vi);

if 0.98*ntt > 1
    umt = mean(ui,4);
    umty = squeeze(mean(umt,2));
    vmt = mean(ui,4);
    vmty = squeeze(mean(vmt,2));
else
    umty = squeeze(mean(ui,2));
    vmty = squeeze(mean(vi,2));
end

ups = 0*ones(ix,iz);
%vps = ups;
uv = ups;
wv = ups;
wu = ups;



for ix = 1:nxrr
        for iz = 1:nz
            up = ui(ix,:,iz,:) - umty(ix,iz);
            vp = vi(ix,:,iz,:) - vmty(ix,iz);
            wp = wi(ix,:,iz,:)  - wmean(ix,iz);
            if 0.98*ntt > 1
                uu = mean(up.*up,4);
                ups(ix,iz) = squeeze(mean(uu,2));
                %uu = mean(vp.*vp,4);
               % vps = mean(uu,2);
                uu = mean(up.*vp,4);
                uv = squeeze(mean(uu,2));
                uu = mean(wp.*up,4);
                wu = squeeze(mean(uu,2));
                uu = mean(wp.*vp,4);
                wv = squeeze(mean(uu,2));
            else
                ups(ix,iz) = squeeze(mean(up.*up,2));
               % vps = mean(vp.*vp,2);
                uv(ix,iz) = squeeze(mean(up.*vp,2));
                wu(ix,iz) = squeeze(mean(up.*wp,2));
                wv(ix,iz) = squeeze(mean(vp.*wp,2));
            end
        end
end
                
clear ui
clear vi
clear wi

wuz = 0*ones(nxrr,nz);
wvz = wuz;
uux = 0*wuz;
uvx = uux;
uus = uux;
uvs = uux;



wuz(:,2:end-1) = 0.5*(wu(:,3:end) - wu(:,1:end-2))./dz(:,2:end-1);
wvz(:,2:end-1) = 0.5*(wv(:,3:end) - wv(:,1:end-2))./dz(:,2:end-1);
clear wu
clear wv


cc = sum(umty.*dz.*wuz,2)./hh;
cs = mean(cc);
cc = sum(vmty.*dz.*wvz,2)./hh;
cs = cs + mean(cc);



uus(:,2:end-1) = 0.5*(ups(:,3:end) - ups(:,1:end-2))/ds;
uvs(:,2:end-1) = 0.5*(uv(:,3:end) - uv(:,1:end-2))/ds;
uux(2:end-1,:) = 0.5*(ups(3:end,:) - ups(1:end-2,:))/dx;
uvx(2:end-1,:) = 0.5*(uv(3:end,:) - uv(1:end-2,:))/dx;


uux = uux + uus.*dxc;
uvx = uvx + uvs.*dxc;

cc = sum((vmty.*uvx + umty.*uux).*dz,2)./hh;
cb = mean(cc);




netcdf.close(ncid)