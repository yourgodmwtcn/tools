% TTP_plot.m  2/26/2009  Parker MacCready
%
% plots the output of PS_tides for a section off of Three Tree Point

clear
load('SEG481_090301.mat');

figure

dlo = datenum(2009,3,8,0,0,0);
dhi = datenum(2009,3,14,0,0,0);

% Surface Height
subplot(211)
plot(Time_datenum_PST,Height,'-k');
aa = axis;
axis([dlo dhi aa(3:4)]);
datetick('x',6,'keeplimits')
ylabel('Height (m)');
title('I think that zero = MLLW')
grid on

%Current
subplot(212)
plot(Time_datenum_PST,Current,'-k');
aa = axis;
axis([dlo dhi aa(3:4)]);
datetick('x',6,'keeplimits')
xlabel(['Time ',num2str(jday_year),' PST']);
ylabel('U (m s^{-1})');
title('Section-Averaged Current (FLOOD POSITIVE)');
grid on
