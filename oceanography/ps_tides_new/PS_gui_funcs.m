function PS_gui_funcs(type, lookup)

global PS_AREA CHAN SEG_INDEX lat_sh lon_sh TIDE_OUT
global h0 hAbx hCbx hSbx hSG hDT hDY hMap hCalc

switch type
    case 1
        ARind = get(hAbx, 'Value');
        CHnm = PS_AREA(ARind).chan_name;
        set(hCbx, 'Value',1, 'String',CHnm);
        set(hSbx, 'Value',1, 'String',[]);
        set(hSG, 'String',[], 'UserData',NaN);
        PS_gui_plot(1, 0);
    case 2
        CHind = get(hCbx, 'Value'); ARind = get(hAbx, 'Value');
        xx = PS_AREA(ARind).chan_no;
        yy = CHAN(xx(CHind)).segment;
        set(hSbx, 'Value',1, 'String',num2str(yy));
        set(hSG, 'String',[], 'UserData',NaN);
        PS_gui_plot(1, 0);
    case 3
        if lookup>0
            SEGent = [get(hSG, 'String') ' ']; % make sure not empty
            SEGent = str2num(SEGent); isn = [];
            if ~isempty(SEGent)
                isn = find(SEG_INDEX.segno == SEGent);
            end
            if ~isempty(isn)
                ARind=SEG_INDEX.Area(isn);
                CHnum=SEG_INDEX.Channo(isn);
                CHnm=PS_AREA(ARind).chan_name;
                CHno=PS_AREA(ARind).chan_no;
                CHv=find(CHno==CHnum);
                xx=CHAN(CHnum).segment; SEGv=find(xx==SEGent);
                % update list_boxes
                set(hAbx, 'Value',ARind);
                set(hCbx, 'Value',CHv, 'String',CHnm, 'ListboxTop',max([CHv-1 1]));
                set(hSbx, 'Value',SEGv, 'String',num2str(xx), ...
                    'ListboxTop',max([SEGv-1 1]));
                set(hSG, 'UserData',1);
                PS_gui_plot(3-lookup, 1);
            else
                set(hSG, 'String', [num2str(SEGent) ' not found'], 'UserData',NaN);
            end
        else
            SEGind = get(hSbx, 'Value'); xx = str2num(get(hSbx, 'String'));
            SEGent = xx(SEGind);
            set(hSG, 'String',num2str(SEGent), 'UserData',1);
            PS_gui_plot(2, 1);
        end
    case 4
        axes(hMap); xy = get(hMap, 'CurrentPoint'); xy=xy(1,1:2);
        xl = get(hMap, 'XLim'); yl = get(hMap, 'YLim');
        xs = SEG_INDEX.lonpts(:,2); ys = SEG_INDEX.latpts(:,2); % segment midpoints
        iSgs = find(xs>=xl(1) & xs<=xl(2) & ys>=yl(1) & ys<=yl(2)); % segs in plot
        dx = xs(iSgs)-xy(1); dy = ys(iSgs)-xy(2);
        dxy = dx.^2+dy.^2; % distance^2 to clicked point
        [xx ix] = min(dxy); % closest segment
        isn = iSgs(ix); SEGent = SEG_INDEX.segno(isn);
        set(hSG, 'String',num2str(SEGent));
        % Use nearest segment as if entered in text box
        PS_gui_funcs(3, 2-lookup); % 1=zoom to channel, 2=zoom to area
    case 5
        DTstr = get(hDT, 'String'); dterr=0;
        if length(DTstr>4)
            eval(['dtn = datenum(DTstr);'], 'dterr=1;');
            if ~dterr
                dts = datestr(dtn,1); % re-format to verify date
                set(hDT, 'String',dts, 'UserData',1);
            end
        else
            dterr=1;
        end
        if dterr
            set(hDT, 'String','Enter VALID date', 'UserData',NaN);
        end
    case 6
        DYnum = [get(hDY, 'String') ' ']; % make sure not empty
        DYnum = str2num(DYnum); DYint = round(DYnum);
        if isempty(DYnum) | DYint~=DYnum | DYint<1
            set(hDY, 'String','1 to 450', 'UserData',NaN);
        elseif DYint>450
            set(hDY, 'String','1 to 450', 'UserData',NaN);
        else
            set(hDY, 'String',num2str(DYint), 'UserData',1);
        end
    case 7
        SEGent = get(hSG, 'String'); SEGok = get(hSG, 'UserData');
        DTent = get(hDT, 'String'); DTok = get(hDT, 'UserData');
        DYent = get(hDY, 'String'); DYok = get(hDY, 'UserData');
        if isnan(SEGok + DTok + DYok)
            errordlg('Segment, Start Date, and No. of Days are REQUIRED', ...
                'Puget Sound Tidal Computation');
        else
            PS_gui_compute( str2num(SEGent), DTent, str2num(DYent) );
        end

end  % of switch type

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PS_gui_plot(pzoom, segHI)

global PS_AREA CHAN SEG_INDEX lat_sh lon_sh
global h0 hAbx hCbx hSbx hSG hDT hDY hMap hCalc

ARind = get(hAbx, 'Value');
xx = PS_AREA(ARind).chan_no;
CHind = get(hCbx, 'Value');
SEGch = CHAN(xx(CHind)).segment;
SEGent = [get(hSG, 'String') ' ']; % make sure not empty
SEGent = str2num(SEGent);
iSx = SEG_INDEX.segind; % indices by seg#
Slon = SEG_INDEX.lonpts; Slat = SEG_INDEX.latpts;

switch pzoom
    case 1
        lims = SEG_INDEX.LonLatLims(ARind,1:4) + [-1 1 -1 1]/20;
    case 2
        xx = Slon(iSx(SEGch),:); yy = Slat(iSx(SEGch),:);
        % include extent of channel
        lims = [min(min(xx)) max(max(xx)) min(min(yy)) max(max(yy))] ...
            + [-1 1 -1 1]/50;
        % make at least 0.1-by-0.1 degrees
        dx = 0.1 - (lims(2)-lims(1)); dy = 0.1 - (lims(4)-lims(3));
        if dx>0
            lims([1 2]) = lims([1 2]) + [-dx/2 dx/2];
        end
        if dy>0
            lims([3 4]) = lims([3 4]) + [-dy/2 dy/2];
        end
end

axes(hMap)
cla
plot(lon_sh, lat_sh, 'k-', 'linewidth',2, ...
    'ButtonDownFcn',['PS_gui_funcs(4,1);']);
axis(lims); hold on
set(hMap, 'DataAspectRatio', [1 cos(mean(lims(3:4))*pi/180) 1]);

% Find indices of segments within plot limits
xs = Slon(:,2); ys = Slat(:,2); % segment midpoints
iSgs = find(xs>=lims(1) & xs<=lims(2) & ys>=lims(3) & ys<=lims(4));
% Plot segments, with color showing:
%   out-of-area, out-of-channel, in-channel, chosen segment
pcol='gbmr';
for ic=1:length(iSgs)
    sn = SEG_INDEX.segno(iSgs(ic));
    ip=1; lth=1; zt='0';
    if SEG_INDEX.Area(iSgs(ic))==ARind
        ip=2;
        if ~isempty( find(SEGch==sn) )
            ip=3; lth=1.3; zt='1';
            if ~isempty(SEGent) & segHI & SEGent==sn
                ip=4; lth=2.5;
            end
        end
    end
    plot(Slon(iSgs(ic),:), Slat(iSgs(ic),:), [pcol(ip) '-'], ...
        'LineWidth',lth, 'ButtonDownFcn',['PS_gui_funcs(4,' zt ');']);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PS_gui_compute( Seg, DTbeg, Ndays )

global PS_AREA CHAN SEG_INDEX lat_sh lon_sh TIDE_OUT
global h0 hAbx hCbx hSbx hSG hDT hDY hMap hCalc

Heightm =  PStide_at_seg( Seg, DTbeg, Ndays, 1 );
Currentms =  PStide_at_seg( Seg, DTbeg, Ndays, 2 );

tmp = datevec(DTbeg);
IYear = tmp(1); % 4-digit year
ya = num2str(10000+IYear); ya=ya(4:5);
ma = num2str(100+tmp(2)); ma=ma(2:3);
da = num2str(100+tmp(3)); da=da(2:3);
filnam = ['SEG' num2str(Seg) '_' ya ma da '.mat'];

elJ85 = datenum('01-Jan-1985');
elJyr = datenum(IYear,1,1);
eldt = datenum(DTbeg);
IYd = eldt - elJyr + 1; % Julian day for year (Jan 1, IYear = 1)
IJd = eldt - elJ85 + 1; % Elapsed days since Jan 1, 1985 (Jan 1, 1985 = 1)

jday_PST = IYd + ([0:length(Heightm)-1]' +0)/24; % +8 if UTS instead
Segment = Seg;
isn = SEG_INDEX.segind(Seg);
LatLon = [SEG_INDEX.latlon1(isn,:) SEG_INDEX.latlon2(isn,:)];
Channel = [num2str(SEG_INDEX.Channo(isn)) ' ' SEG_INDEX.chan{isn}];

nt = length(TIDE_OUT) + 1;
TIDE_OUT(nt).Height = Heightm;
TIDE_OUT(nt).Current = Currentms;
TIDE_OUT(nt).Segment = Segment;
TIDE_OUT(nt).Channel = Channel;
TIDE_OUT(nt).jday_year = IYear;
TIDE_OUT(nt).jday_PST = jday_PST;
TIDE_OUT(nt).LatLon = LatLon;
TIDE_OUT(nt).StartDatePST = DTbeg;
% PM Edit 2/26/2009 Add a datenum format time output field
Tdatenum = datenum(IYear,1,1,0,0,0) + jday_PST - 1;
TIDE_OUT(nt).Time_datenum_PST = Tdatenum;
TIDE_OUT(nt).Time_datenum_GMT = Tdatenum + 8/24;

TIDE = TIDE_OUT(nt);

PSask = 1;
while PSask
    ActStr = questdlg('Select action for this data:', ...
        ['Segment ' num2str(Seg) ' Tidal Data'], 'Plot', 'Save', 'Next', 'Next');
    switch ActStr
        case 'Plot'
            if 0
                figure
                subplot(2,1,1)
                plot(jday_PST, Currentms);
                xlim([jday_PST(1) jday_PST(end)]);
                ylabel('Current (m s^{-1})'); grid on
                title(['Segment ' num2str(Seg) ' start ' DTbeg '(PST) [FLOOD POSITIVE]']);
                subplot(2,1,2)
                plot(jday_PST, Heightm);
                xlim([jday_PST(1) jday_PST(end)]);
                ylabel('height / m'); grid on
                xlabel(['Julian day, ' num2str(IYear) ' (PST)']);
            else
                figure
                % Surface Height
                subplot(211)
                plot(Tdatenum,Heightm,'-k');
                datetick('x',6,'keeplimits')
                ylabel('Height (m)');
                grid on
                %Current
                subplot(212)
                plot(Tdatenum,Currentms,'-k');
                datetick('x',6,'keeplimits')
                xlabel(['Time ',num2str(IYear),' PST']);
                ylabel('U (m s^{-1})');
                title('Section-Averaged Current (FLOOD POSITIVE)');
                grid on
            end
        case 'Save'
            [sfil,spath] = uiputfile(filnam,'Save file name');
            if sfil~=0
                fstr=['save ''' spath sfil ''''];
                vars=fieldnames(TIDE);
                for ic=1:length(vars)
                    eval([vars{ic} '= TIDE.' vars{ic} ';']);
                    fstr = [fstr ' ' vars{ic}];
                end
                eval(fstr)
                PSask=2;
            end
        case 'Next'
            if PSask<2
                helpdlg(['Data are in structure TIDE_OUT(' num2str(nt) ')'], ...
                    'Finished with Tidal Segment');
            end
            PSask = 0;
    end % of switch PSask
end % of while


return
