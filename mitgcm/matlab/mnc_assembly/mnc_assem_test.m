%
%  Ed Hill
%  Fri Nov 26 20:08:14 EST 2004
%

%  Tests of MNC output assembly and disassembly.

%  matlab -nojvm

clear all
close all

disp(sprintf('\nStarting test with "exch1" geometry:\n'));

vars = struct([]);
vars(1).name = 'iter';
vars(2).name = 'U';
vars(3).name = 'Unk';
vars(4).name = 'V';
vars(5).name = 'Temp';
vars(6).name = 'S';

fpat = 'exp0_20041126_0001/state.0000.%06d.nc';
[nt,nf] = mnc_assembly(fpat,vars);

ncload('all.00000.nc')
figure(1)
subplot(2,2,1),  a = squeeze(Temp(2,1,:,:)); surf(a), view(2)
subplot(2,2,2),  a = squeeze(U(2,1,:,:));    surf(a), view(2)
subplot(2,2,3),  a = squeeze(V(2,1,:,:));    surf(a), view(2)
subplot(2,2,4),  a = squeeze(S(2,1,:,:));    surf(a), view(2)

disp(sprintf('\nTest with "exch1" geometry complete.\n'));


disp(sprintf('\nStarting test with "exch2" geometry:\n'));

fpat = 'aim.5l_cs_20041209_0001/state.0000.%06d.nc';
[nt,nf] = mnc_assembly(fpat,vars);

ncload('all.00000.nc')
figure(2)
for i = 1:6
  subplot(6,4,(i-1)*4+1)
  a = squeeze(Temp(2,1,i,:,:)); surf(a), view(2)
  subplot(6,4,(i-1)*4+2)
  a = squeeze(U(2,1,i,:,:));    surf(a), view(2)
  subplot(6,4,(i-1)*4+3)
  a = squeeze(V(2,1,i,:,:));    surf(a), view(2)
  subplot(6,4,(i-1)*4+4)
  a = squeeze(S(2,1,i,:,:));    surf(a), view(2)
end

disp(sprintf('\nTest with "exch2" geometry complete.\n'));
