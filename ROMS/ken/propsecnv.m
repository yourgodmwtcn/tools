
function propsecnv(filename,yfar,tlow,thigh,xp)

% plot sections of v and T at location = yfar using the new netcdf
%   for stability runs
%   Does section of Ertel Vorticity
%   Call as
%       propsecnv(filename,yfar,tlow, thigh,xp)
%   where
%       filename is a string, e.g.'XSHELFCI003_avg.nc''
%       yfar is in km
%       times are in days
%       xp = x position to profile in km

%   KHB 7/9/2008, 3/5/2010  1/4/2011  6/7/2011  8/27/2011
%   correct vertical grid 8/29/11

rhozero = 1027;
tzero = 14.0;
beta = 1.7e-4;



ncid = netcdf.open(filename,'NC_NOWRITE');
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'))/86400;
ntime = length(time);
tchar = filename(9:11);

vp = 1;         % working with outflow (v > 0)
%vp = -1;       % working with inflow (v < 0)

rt = filename(7:8);

tmax = max(time);
if tmax < tlow
    disp(['max model time is less than start time, tmax = ' num2str(tmax)])
    tlow = tmax;
end
itime = find(time <= thigh & time >= tlow);
itiml = min(itime);
ntt = length(itime);

srho = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'s_rho'));
nz = length(srho);
csr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Cs_r'));
theta_s = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_s'));
theta_b = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'theta_b'));
Tcline = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Tcline'));
vxform = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vtransform'));
vstretch = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'Vstretching'));
hc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hc'));

user = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'user'));

xzero = user(5);


xv = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_v'));
[nxv,nyv] = size(xv);
xv = squeeze(xv(:,end))/1000;
%dxv = xv(2)-xv(1)
yv = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y_v'),[10,0],[1,nyv])/1000;
yv = squeeze(yv);

xu = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_u'));
[nxu,nyu] = size(xu);
xu = squeeze(xu(:,end))/1000;
yu = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y_rho'),[10,0],[1,nyu])/1000;
yu = squeeze(yu);

xr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_rho'));
[nxr,nyr] = size(xr);
xr = squeeze(xr(:,end))/1000;
yr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y_rho'),[10,0],[1,nyr])/1000;
yr = squeeze(yr);
dy = yr(2) - yr(1);

mr = nxr;
nr = nyr;
[xtau,itau] = min(abs(xr-xzero));
dx = (xr(2)-xr(1));
ipro = round(xp/dx + 0.5);
iproul = ipro-1;
itau = ipro;

[yfara,ifar] = min(abs(yr-yfar));
yfara = yr(ifar);
h = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'),[0,ifar-1],[nxr,1]);
h = squeeze(h);

hx = NaN*h;
hx(2:end-1) = 0.5*(h(3:end) - h(1:end-2))/dx;
hx(1) = (h(2) - h(1))/dx;
hx(end) = (h(end) - h(end-1))/dx;
f = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'f'),[2,ifar-1],[1,1]);
f = squeeze(f);

% wall mask


    xfw = NaN*[1 1 1];
    zfw = xfw;
    
dy = dy*1000;
dx = dx*1000;

temp = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),[0,ifar-2,0,itiml-1],[nxr,3,nz,ntt]);
temp = squeeze(temp);
tempy = squeeze(0.5*(temp(:,3,:,:)-temp(:,1,:,:))/dy);
temp = squeeze(temp(:,2,:,:));



vread = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0,ifar-1,0,itiml-1],[nxv,2,nz,ntt]);

if length(itime) > 1.5
    v = 0.5*(vread(:,1,:,:) + vread(:,2,:,:));
else
    v = 0.5*(vread(:,1,:) + vread(:,2,:));
end
v = squeeze(v);
clear vread



uread = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[0,ifar-2,0,itiml-1],[nxu,3,nz,ntt]);
u = squeeze(uread);
clear uread
ui = NaN*ones(mr,3,nz,ntt);
ui(2:end-1,:,:,:) = 0.5*(u(1:end-1,:,:,:) + u(2:end,:,:,:));
uy = squeeze(0.5*(ui(:,3,:,:) - ui(:,1,:,:))/dy);
u = squeeze(u(:,2,:,:));
ui = squeeze(ui(:,2,:,:));
%if isempty(imasku) ~= 0
%    u(:,:,imasku) = u(:,:,imasku)*NaN;
%end

vb = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'v'),[0,ifar-1,0,0],[nxu,1,1,ntime]);
vb = squeeze(vb);
rhob = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'rho0'));
%cd = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'rdrg2'));



hm = max(max(h));

tpro = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),[ipro-1,ifar-1,0,0],[1,1,nz,ntime]);
tpro = squeeze(tpro);
upro = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[iproul-1,ifar-1,0,0],[2,1,nz,ntime]);
if ntt > 1.5
uproo = 0.5*(upro(:,:,1,1) + upro(:,:,1,2));
uproo = 0.5*(upro(1,1,:,:) + upro(2,1,:,:));
upro = squeeze(uproo);
clear uproo
zpro = hc*srho + (h(ipro)-hc)*csr;



figure(6)
%pcolor(time,zpro,tpro')
%shading interp
%colorbar
ci = 0.1*(-200:200);
%ci = 0:200;
[C,hh] = contour(time,zpro,tpro,ci,'Color','r');
%clabel(C,hh);
%get(gca)

hold on
[C,hh] = contour(time, zpro,upro);
clabel(C,hh)
xlabel('Time (days)')
ylabel('Depth (m)')
title(['T and u at x, y = ' int2str(xr(ipro)) '  ' int2str(yr(ifar)) ' km'])
axis([0 max(time) -h(ipro) 0]);
hold off
clear upro
clear tpro

if rt == 'zp'
%if rt ~= '2D'
    figure(7)
    itm = max(itime);
    tysec = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temp'),[ipro-1,0,0,itm-1],[1,nyr,nz,1]);
    tysec = squeeze(tysec);
    tysec(1,:) = NaN;
    tysec(end,:) = NaN;
    uysec = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'u'),[ipro-1,0,0,ntime-1],[1,nyu,nz,1]);
    uysec = squeeze(uysec);

      
    [C,hh] = contour(yr,zpro,tysec',ci,'Color','r');
    clabel(C,hh)

    %hold on
    %[C,hh] = contour(yu,zpro,uysec,'Color','b');
    %clabel(C,hh);
    %hold off
    xlabel('y (km)')
    ylabel('z (m)')
    title(['Temperature at x, t = ' int2str(xr(ipro)) ',  ' int2str(time(itm))])
    clear tysec
    clear uysec
    %axis([0 yr(end) -h(ipro) (-h(ipro)+20)])
end
    

end



iyts = round(nr*0.5);
ixts = round(mr*0.8);
iyts = ifar;
ixts = ipro;


vts = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vbar'),[ixts-1,iyts-1,0],[1,1,ntime]);
tauy = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bvstr'),[ixts-1,iyts-1,0],[1,1,ntime]);
tauy = squeeze(tauy);


figure(3)
subplot(2,1,1)
plot(time,squeeze(vts))
xlabel('Time (days)')
ylabel('vbar (m/sec)')
title(['vbar at x,y = ' int2str(xr(ixts)) ',  ' int2str(yr(iyts)) ])
subplot(2,1,2)
plot(time,tauy,[0 max(time)],[0 0],'k')
xlabel('Time(days)')
ylabel('y stress (N/m^2)')
title(['y stress at x,y, = ' int2str(xr(itau)) ',  ' int2str(yfar)])


xrf = NaN*ones(mr,nz);
zrf = xrf;
dxc = xrf;
ds = srho(2) - srho(1);
dz = xrf;
zs = xrf;

for jj = 1:mr
  
    
    xrf(jj,:) = xr(jj)*ones(nz,1);
    if round(vxform) ==1
        zrf(jj,:) =  hc*srho + (h(jj)-hc)*csr;
    else
        sss = (hc*srho + h(jj)*csr)/(hc+h(jj));
        zrf(jj,:) = h(jj)*sss;;
    end
    
    zs(jj,2:end-1) = 0.5*(zrf(jj,3:end) - zrf(jj,1:end-2))/ds; 
    zs(jj,1) = (0.5*(zrf(jj,2)+zrf(jj,1)) +h(jj))/ds;
    zs(jj,end) = -0.5*(zrf(jj,end-1)+zrf(jj,end))/ds;
    dxc(jj,:) = -hx(jj)*csr./zs(jj,:)';
    dz(jj,2:end-1) = (zrf(jj,3:end)-zrf(jj,1:end-2))/2;
    dz(jj,1) = (zrf(jj,2)+zrf(jj,1))/2 + h(jj);
    dz(jj,end) = -(zrf(jj,end) + zrf(jj,end-1))/2;
    
    
end

%figure(11)
%plot(dz(10,:),'+')
%figure(12)
%plot(zrf(10,:),'+')



%xuf = 0*squeeze(u(1,:,:));
%zuf = xuf;
%size(h)
xuf = 0*ones((mr-1),nz);
zuf = xuf;
[nux,nuz] = size(xuf);




for jj = 1:nux
    xuf(jj,:) = xu(jj)*ones(nuz,1);
     if round(vxform) ==1
        zuf(jj,:) =  hc*srho + ((h(jj)*h(jj+1))*0.5 - hc)*csr;
    else
        sss = (hc*srho + h(jj)*csr)/(hc+0.5*(h(jj)+h(jj+1)));
        zuf(jj,:) = h(jj)*sss;;
    end
    %zuf(jj,:) = hc*srho + ((h(jj)+h(jj+1))*0.5 - hc)*csr;  old form
end

%   Work out bottom mask
xf = [xr' max(xr) min(xr)];
hmaxx = max(h);
zf = [-h(:)' -hmaxx -hmaxx];

hm = max(h);
  if round(vxform) ==1
        zu =  hc*srho + (hm - hc)*csr;
    else
        sss = (hc*srho + hm*csr)/(hc+hm);
        zu = h(jj)*sss;;
    end
%zu = hc*srho + (hm-hc)*csr;    old form

xlow = 0.5*(xr(1) + xr(2));
xhigh = 0.5*(xr(end) + xr(end-1));
axx = [xlow xhigh -hmaxx 0];



ntt = length(itime);


for iii = 1:ntt
    ii = itime(iii);
    tt = time(ii);
    %pcolor(xrf,zrf,squeeze(temp(iii,:,:)))
    %shading interp
    %colorbar
    %caxis([10 20])
    
    sqt = squeeze(temp(:,:,iii));
    sqty = squeeze(tempy(:,:,iii));
    sqv = squeeze(v(:,:,iii));
    squ = squeeze(u(:,:,iii));
    squi = squeeze(ui(:,:,iii));
    squy = squeeze(uy(:,:,iii));
   if length(itime) < 1.5
       sqv = v;
       sqt = temp;
       squ = u;
       squi = ui;
       sqty = tempy;
       squy = uy;
   end
   
   
    iim = find(abs(sqv) > 2.0);
    if isempty(iim) == 0
        sqv(iim) = NaN;
    end
    vmaxmax = max(max(sqv));
    vminmin = min(min(sqv));
    vcontmax = 0.1*vmaxmax;
    vcont = 0.1*vminmin;
    vcont20 = 0.2*vminmin;
    vmcont20 = 0.2*vmaxmax;
    disp(['Extreme v at t = ' num2str(tt) ' are ' num2str(vminmin) ',  ' num2str(vmaxmax) ' m/sec'])
    
    cit = -30:30;
    
    
    figure(1)
    contour(xrf,zrf,sqt,cit,'r')
    ylabel('Depth (m)')
    xlabel('x (km)')
    title(['Temperature and v at y, t = ' num2str(yfar) ' km, ' num2str(tt) ' days'])
    
    hold on
    [C,hhh] = contour(xrf,zrf,sqv,'k');
    clabel(C,hhh);
  
    fill(xf,zf,0.9*[1 1 1])
   
    hold off
    axis(axx)

  
    
    
    %       Plot surface velocity vs. x
    figure(10)
    plot(xrf(:,end),sqv(:,end))
    hold on
    plot(xrf(:,end),sqv(:,end)*0, 'k')
    hold off
    title('Surface velocity')
    xlabel('x (km)')
    
    figure(4)
   % [C,hh] = contour(xuf,zuf,squeeze(u(iii,:,:)),'k');
   % clabel(C,hh)
   
  
    iim = find(abs(squ) > 2.0);
    if isempty(iim) == 0
        squ(iim) = NaN;
    end
    pcolor(xuf,zuf,double(squ))
    shading interp
     colorbar
    
    xlabel('x (km)')
    ylabel('z (m)')
    title(['u (m/sec) at y = ' int2str(yu(ifar)/1000) ' km, ' num2str(tt) ' days'])
    hold on
    fill(xf,zf,0.9*[1 1 1])
   %fill(xfw,zfw,0.9*[1 1 1])
    hold off
    axis(axx)
    
    
   
  %%    Plot Ertel Vorticity  
   figure(2)
    rho = -rhozero*beta*(sqt-tzero);
    rhoz = NaN*rho;
    rhox = NaN*rho;
    rhoy = -rhozero*beta*sqty;
    %   I have uy already
    vx = rhox;
    vz = rhox;
    uz = rhox;
    
    dx = (xrf(3,1) -xrf(2,1))*1000;
    
    vz(:,2:end-1) = 0.5*(sqv(:,3:end) - sqv(:,1:end-2))/ds;     % actually, dv/ds
    vz(:,end) = 0;
    rhoz(:,2:end-1) = 0.5*(rho(:,3:end)- rho(:,1:end-2))/ds;
    rhoz(:,end) = NaN;
   % vz is OK
    for jx = 2:mr-1
        vx(jx,:) = 0.5*(sqv(jx+1,:) - sqv(jx-1,:))/dx;
        vx(jx,:) = vx(jx,:) +  dxc(jx,:).*vz(jx,:);
        rhox(jx,:) = 0.5*(rho(jx+1,:) - rho(jx-1,:))/dx;
        rhox(jx,:) = rhox(jx,:) +  dxc(jx,:).*rhoz(jx,:);
    end
   
   
    for jx = 1:mr
       
        uz(jx,2:end-1) = 0.5*(squi(jx,3:end) - squi(jx,1:end-2))./dz(jx,2:end-1);
        uz(jx,end) = 0;
        vz(jx,:) = vz(jx,:)./zs(jx,:);
        rhoz(jx,:) = rhoz(jx,:)./zs(jx,:);
    end
    %
    
    
    ert = -((rhoz.*(f + vx-squy) - rhox.*vz + rhoy.*uz)/rhozero);
   % ert = -rhoz*f/rhozero;
    ert = ert*rhozero/f;
    figure(2)
    dpv = 4e-4;
    cim = -dpv*(0:40)/10;
    [C,hh]=contour(xrf,zrf,ert,cim,'r');
    clabel(C,hh)
    %contour(xrf,zrf,ert,'r');
    cip = 2*dpv*(1:40);
    hold on
        [C,hh]= contour(xrf,zrf,ert,cip,'k');
        clabel(C,hh)
        fill(xf,zf,0.9*[1 1 1])
    hold off
    
    ylabel('Depth (m)')
    xlabel('x (km)')
    title(['Ertel vorticity (times \rho_0/f) at y, t = ' num2str(yfar) ' km, ' num2str(tt) ' days. Note different CI for +/-'])
    axis(axx)

   
%%%%%


   figure(5)
   nsxx = 4;
   temps = 5;
   plot([0 0],[-10 0])
   hold on
   dx = xrf(2,1)-xrf(1,1);
   for iix = 2:nsxx:mr
        tt = sqt(iix,:);
        setoff = dx*(iix-1) - dx/2;
        ttt = setoff + temps*squeeze(tt);
        ttt = ttt - temps*tt(1);
        plot(ttt,zrf(iix,:),'k')
   end
    fill(xf,zf,0.9*[1 1 1])
     [C,hh] = contour(xrf,zrf,sqv,'b');
    clabel(C,hh);
    
     xmax = mr*dx - dx/2;
    zmin = min(min(zrf));
    zmin = -hmaxx;
    
    axis(axx)
    
    zlow = zmin*0.8;
    xlow = xmax*0.1;
    dxt = 2*temps;
    plot([xlow (xlow+dxt)],[zlow,zlow],'k')
    text(xlow,zlow + 10,'2 degrees')
      
    hold off
    title('Temperature profiles and v')
    ylabel('z (m)')
   
   pause
   
end

netcdf.close(ncid)
