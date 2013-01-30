% Calculate PV on C-grid. Result at interior W points
%        [pv,xpv,ypv,zpv] = pv_cgrid(rgrid,u,v,rho,f,rho0)
% supply rgrid structure with 
%       xu,yu,zu & xv,yv,zv & xr,yr,zr (all matrices) 
%                   & zw & s_w(vectors)

function [pv,xpv,ypv,zpv] = pv_cgrid(rgrid,u,v,rho,f,rho0)

    xpv = rgrid.xr(2:end-1,2:end-1,end);
    ypv = rgrid.yr(2:end-1,2:end-1,end);
    zpv = avg1(rgrid.zr(2:end-1,2:end-1,:),3);
    
    gridu.xmat = rgrid.xu; gridu.ymat = rgrid.yu; gridu.zmat = rgrid.zu;
    gridv.xmat = rgrid.xv; gridv.ymat = rgrid.yv; gridv.zmat = rgrid.zv;
    gridr.xmat = rgrid.xr; gridr.ymat = rgrid.yr; gridr.zmat = rgrid.zr;
    
    gridu.s = rgrid.s_rho; gridv.s = rgrid.s_rho; gridr.s = rgrid.s_rho;
    
    s = size(rho);
    if length(s) == 3
        s(4) = 1;
    end
    
    if length(f) == 1 % constant f0
        f   = repmat(f,[1 1 s(3) s(4)]);
    end
    
    % average to land on vx - uy points
    f = avg1(avg1(avg1(f,2),2),1);
        
    if ~exist('rho0','var');
        fprintf('\n\n Using rho0 = 1027. \n\n');
        rho0 = 1027;
    end
    
    % calculate gradients
    vx    = horgrad_cgrid(rgrid,gridv,v,1); 
    vy    = horgrad_cgrid(rgrid,gridv,v,2); 
    vz    = bsxfun(@rdivide,diff(v,1,3),diff(rgrid.zv,1,3));

    ux    = horgrad_cgrid(rgrid,gridu,u,1); 
    uy    = horgrad_cgrid(rgrid,gridu,u,2); 
    uz    = bsxfun(@rdivide,diff(u,1,3),diff(rgrid.zu,1,3));

    tx    = horgrad_cgrid(rgrid,gridr,rho,1); 
    ty    = horgrad_cgrid(rgrid,gridr,rho,2); 
    tz    = bsxfun(@rdivide,diff(rho,1,3),diff(rgrid.zr,1,3));

    % PV calculated at interior rho points
                                % f + vx - uy                      (rho)_z
    pv = -1* double((avgx(avgz(bsxfun(@plus,avg1(vx - uy,2),f)))  .*  tz(2:end-1,2:end-1,:,:) ...
                   - avgy(vz(2:end-1,:,:,:)).*avgz(avgx(tx(:,2:end-1,:,:))) ... % vz * (rho)_x
                   + avgx(uz(:,2:end-1,:,:)).*avgz(avgy(ty(2:end-1,:,:,:))))./rho0);%avgz(lambda(2:end-1,2:end-1,:,:))); % uz*(rho)_y
               
    debug = 0;
    if debug
        pv1 = -avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  tz(2:end-1,2:end-1,:,:);
        pv2 = avgy(vz(2:end-1,:,:,:)).*avgz(avgx(tx(:,2:end-1,:,:)));
        pv3 = avgx(uz(:,2:end-1,:,:)).*avgz(avgy(ty(2:end-1,:,:,:)));
        
        tind = 1;
        yind = 3;
        
        figure;
        contourf(xpv,zpv,squeeze(pv1(:,yind,:,tind))',20);colorbar
        title('(f + v_x -u_y)\rho_z');
        figure;
        contourf(xpv,zpv,squeeze(pv2(:,yind,:,tind))',20);colorbar
        title('v_z \rho_x');
        figure;
        contourf(xpv,zpv,squeeze(pv3(:,yind,:,tind))',20);colorbar
        title('u_x \rho_y');
        figure;
        contourf(xpv,zpv,squeeze(pv(:,yind,:,tind))',20);colorbar
        title('Full PV');
        colormap(hsv);
        pause;
    end
    

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;
    
% following for *flat bottom* only
%     vx    = bsxfun(@rdivide,diff(v,1,1),diff(grid.xv)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     vy    = bsxfun(@rdivide,diff(v,1,2),diff(grid.yv')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
%     ux    = bsxfun(@rdivide,diff(u,1,1),diff(grid.xu)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     uy    = bsxfun(@rdivide,diff(u,1,2),diff(grid.yu')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
%     tx    = bsxfun(@rdivide,diff(rho,1,1),diff(grid.xr)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     ty    = bsxfun(@rdivide,diff(rho,1,2),diff(grid.yr')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    
