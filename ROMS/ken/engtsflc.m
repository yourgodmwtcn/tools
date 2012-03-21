function [tim,eke] = engtsflc(fnameb,nf,dt,runt,xlow,xhigh,npl)

%   Compute time series of energy, heat flux, and PE-KE conversion
%   call as
%       [tim,eke] = engtsflc(fnameb,nf,dt,runt,xlow,xhigh,npl)
%   where%
%       fnameb = basic file name (string:first file to read)
%       nf = number of output files
%       dt = time step (days)
%       runt = total run time (days)
%       xlow,xhigh (km) = x range for computation
%       npl = # points to use for each growth rate estimate

%   6/16/11   8/17/2011  9/22/11

nts = 1            % stride for reading

stt = fnameb(10:12);
if stt == 'avg'
    nf = 1;
    fname1 = fnameb;
elseif nf == 1
    fname1 = fnameb;
else
    if fnameb(10:12) == 'his'
        fname1 = [fnameb(1:12) '_0001.nc'];
    else
        fname1 = [fnameb(1:13) '_0001.nc'];
    end
end

nread = floor(runt/dt);     % total number of averaged blocks

fname1 = fnameb;


ncid = netcdf.open(fname1,'NC_NOWRITE');
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ocean_time'))/86400;
h = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'h'));
h = h(:,5);
zzero = -max(h)/2;
clear h

tmax = max(time);           %  Amount of time in one file 
dtr = (time(2) - time(1))*nts;
nbin = round(dt/dtr)       % # of points/averaged block

ntimes = floor(tmax/dt);    % # of time averaged bins (in first file)

xr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x_rho'))/1000.;
[nxr,nyr] = size(xr);
xr = xr(:,5);
%[xx,irho] = min(abs(xr-xc));
dx =(xr(2) - xr(1));
ixrl = round(1.5 + xlow/dx);
ixrh = round(1.5 + xhigh/dx);
nxrr = length(ixrl:ixrh);
icent = round(mean(ixrl:ixrh));
disp(['Location for flux calculation = ' num2str(xr(icent)) ' km'])
netcdf.close(ncid)



eke= NaN*ones(nread,1);
timm = eke;
tke = eke;
mpe = eke;
tpe = eke;
pep = eke;
zpe = eke;
erf = eke;
mrf = eke;
encontot = eke;
enconm = eke;
cb = eke;
cs = eke;


rbbb = NaN;

fname = fname1;

figure(1)
clf

for itime = 1:nread
   
    tlow = (itime-1)*dt;
    thigh = tlow + dt;
    nfile = 1 + floor((thigh-0.01)/tmax);
    ffname = fnameb(1:12);
%     if fnameb(10:12) ~= 'his'
%         ffname = fnameb(1:13);
%     end
%     
%     if nfile < 9.5 & nf > 1.5
%         fname = [ffname '_000' int2str(nfile) '.nc'];
%     elseif nfile > 9.5 & nf > 1.5
%         fname = [ffname '_00' int2str(nfile) '.nc'];
%     end
    
    
   % timm(itime) = (tlow+thigh)/2;
    fname = fnameb;
    
    [eke(itime),tke(itime),tpe(itime),pep(itime),zpe(itime),ee,rr,ntt,rbbb,tavg]  = engflcal(fname,tlow,thigh,xlow,xhigh,zzero,rbbb,nts);
    erf(itime) = ee;
    mrf(itime) = rr;
    [encontot(itime),enconm(itime),cb(itime),cs(itime)] = enconcal(fname,tlow,thigh,xlow,xhigh,rbbb,nts);
    
    
    
    if ntt < nbin*0.98
        eke(itime) = NaN;
        encontot(itime) = NaN;
        cs(itime) = NaN;
        cb(itime) = NaN;
    end
    timm(itime) = tavg;
    
%     figure(1)
%     plot(timm,eke)
%     xlabel('Time (days)')
%     ylabel('Eddy KE (m /sec)^2')
%     title(['Eddy kinetic energy for ' ffname ' at x = ' int2str(xlow) ' - ' int2str(xhigh) ' km'])

end;

%tpe = tpe*NaN;
mke = tke-eke;
teng = tpe+tke+zpe;
%epe = tpe-mpe;

plot(timm,eke,'b',timm,mke,'b--')
hold on
    xx = [0 runt];
    plot(xx,0*xx,'k')
    plot(timm,pep,'r',timm,tpe,'r--')
    plot(timm,teng,'k--')
hold off
xlabel('Time (days)')
ylabel('Energy (m /sec)^2')
title(['Energy for ' ffname ' at x = ' int2str(xlow) ' - ' int2str(xhigh) ' km. Blue = KE, Red = PE'])

figure(2) 
plot(timm,eke)

figure(3)
plot(timm,erf,'r',timm,mrf,'b--',timm,0*erf,'k')
xlabel('Time (days)')
ylabel('Averaged density flux')
title(['Density flux for ' ffname ' at x = ' int2str(xlow) ' - ' int2str(xhigh) ' km. Blue = mean, Red = total'])
meanmean = mean(mrf);
meantot = mean(erf);
disp([' Mean total and mean fluxes = ' num2str(meantot) ',  ' num2str(meanmean)])


figure(5)
plot(timm,encontot,'r',timm,enconm,'b', timm,cb,'g',timm, cs,'g--',timm,0*timm,'k')
xlabel('Time (days)')
ylabel('Energy conversion (m^2/s^3)')
title('Energy conversions (Positive for PE to KE or mean to eddy)')
%text(timm(end)*0.8, 0,'Green = barotropic conversion, green dashed = shear conversion')
legend('Total PE to KE','Mean PE to KE','Barotropic MKE to EKE','Shear MKE to EKE', ' ','Location','NorthWest')

%%%%%%%%%%
% Now look at growth rates

[xx,jmax] = max(eke);



tfit = NaN*timm;
tsfit = NaN*timm;
ett = tsfit;
ett(2:end-1) = eke(1:end-2) - 2*eke(2:end-1) + eke(3:end);

efit = tsfit;
afit = tsfit;
bfit = tsfit;
nfits = length(1:jmax);
ekelog = log(eke);


iin = find(ett < 0);
%ekelog(iin) = NaN;

for ii = 1:nfits
    jlow = ii;
    jhigh = ii+npl-1 ;
    
    if jhigh <= nfits
    
        AA = 0*ones(2,2);
        BB = 0*ones(2,1);
    
        AA(1,1) = sum(timm(jlow:jhigh).*timm(jlow:jhigh));
        AA(2,1) = sum(timm(jlow:jhigh));
        AA(1,2) = AA(2,1);
        AA(2,2) = npl;
        BB(1) = sum(ekelog(jlow:jhigh).*timm(jlow:jhigh));
        BB(2) = sum(ekelog(jlow:jhigh));
    
        CC = AA\BB;
        tfit(ii) = mean(timm(jlow:jhigh));
    
        A = CC(1);
         B = CC(2);
         if A < 0;
              A = NaN;
         end
    
    
        fcheck = A*timm(jlow:jhigh) + B;
         err = sum((ekelog(jlow:jhigh) - fcheck).^2);
         efit(ii) = err/npl;
         afit(ii) = A;
         bfit(ii) = B;
    
    
    end
end



[xx,jj] = min(efit);
[xx,jj] = max(afit);
jlow = jj;
jhigh = jj+npl-1 ;
fcheck = afit(jj)*timm(jlow:jhigh) + bfit(jj);
tcheck = timm(jlow:jhigh);
fcheck = exp(fcheck);

figure(2)
hold on 

plot(tcheck,fcheck,'g')

%afit
hold off

tim = timm;

%bfit(jj)
%exp(bfit(jj))

disp(['Largest rate (1/day) = ', num2str(afit(jj)) , ' at time ' num2str(tfit(jj)) ' days'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

EKE = eke;
PE = tpe;
t_en = timm*86400;
MKE = mke;
A = afit;
time_A = tfit*86400;

save('energy-ken.mat','EKE','PE','MKE','t_en','A','time_A');
