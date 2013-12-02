function [gval,labs]=axdate(varargin);
% AXDATE Does nice date labelling of current axis
%   AXDATE with no arguments auto-labels the axis.
%   AXDATE(APPROXTICKS) puts approximately APPROXTICKS
%   labels. Ticks are chosen "nicely" (i.e. in various
%   units of exact hours/days/months) so that the
%   exact number may be a little different. If not
%   specified the existing number of ticks is used
%   as a default for APPROXTICKS.
%
%   AXDATE can also generate an approximate number of
%   ticks for a given time interval and return the label
%   and tick locations without plotting them using:
%   [LABS,LOCS]=AXDATE(APPROXTICKS,TLIMS), where
%   TLIMS=[T_MIN T_MAX]. This is useful when annotation
%   must be produced for some other purpose.
%  
%   AXDATE('y',...) changes the y axis. AXDATE('x',...
%   can be used to specify the x axis explicitly.

% R. Pawlowicz 6/Sep/99
  
ax='x';
if length(varargin)>0 & isstr(varargin{1}),
    ax=lower(varargin{1}(1));
    varargin(1)=[];
end;
    
  
if length(varargin)==2,
  nticks=varargin{1};
  dlim=varargin{2};
elseif length(varargin)==1,
  nticks=varargin{1};
  dlim=get(gca,[ax 'lim']);       % Get curent axis setting
else
  nticks=length(get(gca,[ax 'tick']));
  dlim=get(gca,[ax 'lim']);       % Get curent axis setting
end;
gmin=dlim(1);gmax=dlim(2);


exactint=diff(dlim)/(nticks-1); % Exact interval in days
                                % for a given number of ticks

% These are the intervals which we will allow (they are "nice" 
% in the sense that they come to various nice intervals of days

niceints=[ (1/86400)*[5 6 10 12 15 20 30] ...
           (1/1440)*[1 2 3 4 5 6 10 12 15 20 30] ...
           (1/24)*[1 2 3 4 6 8 12] ...
	   1 2 3 4 5 7 10 14 ...
           (365.25/12)*[1 2 3 4 6] ...
	   365.25*[1 2 3 4 5 10 15 20]];

% If intervals are days are more frequent than everything is nice,
% but when they are months or years we have to fiddle the spacing
% more accurately because of the varying lengths - this is a flag
% indicator for when that is necessary.

intflag=[           0 0 0 0 0 0 0 ... 
                    0 0 0 0 0 0 0  0  0  0  0 ...
                   0 0 0 0 0 0 0 ...
	   0 0 0 0 0 0 0  0 ...
	                1 1 1 1 1 ...
	           2 2 2 2 2 2  2  2];
		   	   
[dun,I]=min(abs(niceints-exactint));


% Calculate the exact tick locations. We must
% shift intervals a little when dealing with month or year
% increments (to start things on an even year and to line them
% up with the first day of months).

switch intflag(I),
  case 0,
    gval=niceints(I)*[ceil(gmin/niceints(I)):fix(gmax/niceints(I))];
 
  case 1,      % month-based ticks
    [Y,M,D,H,MI,S]=datevec(gmin);
    gval=[datenum(Y,1,1):niceints(I):gmax];
    [Y,M,D,H,MI,S]=datevec(gval);
    M=round(M+D/31);
    gval=datenum(Y,M,1);    

  case 2,      % year-based ticks.
    [Y,M,D,H,MI,S]=datevec(gmin);
    gval=[datenum(Y,1,1):niceints(I):gmax];
    [Y,M,D,H,MI,S]=datevec(gval);
    Y=round(Y+M/12);
    gval=datenum(Y,1,1);
end;

% Get rid of labels outside axis

gval=gval(gval>=gmin & gval<=gmax);


% Now we have to pick the labels. Different formats are used
% depending on the length of the axis, but in general I
% try to keep the labels as small as possible, but, e.g. the
% first day of a new month or year is labelled with that
% new day or year, and one of the labels has the year as
% well. If your tastes differ...modify away!

if niceints(I)<1/1440,  % Less than 1 minute
   labs=cellstr(datestr(gval+1/86400,'HH:MM:SS'));
   kk=find(rem(gval,1)==0);    % Things at exact days
   if any(kk),
     labs(kk)=cellstr(datestr(gval(kk),'mm/dd'));
     onelab=datestr(gval(kk(1)),'mm/dd/yy');
     labs{kk(1)}=onelab([7 8 3 1 2 3 4 5]);
   else
     onelab=datestr(gval(1),'mm/dd/yy');
     labs{1}=[onelab([7 8 3 1 2 3 4 5]) ' ' datestr(gval(1),'HH:MM:SS')];
   end;

elseif niceints(I)<1,              % Ticks less than 1 day apart
   labs=cellstr(datestr(gval+1/86400,'HH:MM'));
   kk=find(rem(gval,1)==0);    % Things at exact days
   if any(kk),
     labs(kk)=cellstr(datestr(gval(kk),'mm/dd'));
     onelab=datestr(gval(kk(1)),'mm/dd/yy');
     labs{kk(1)}=onelab([7 8 3 1 2 3 4 5]);
   else
     onelab=datestr(gval(1),'mm/dd/yy');
     labs{1}=[onelab([7 8 3 1 2 3 4 5]) ' ' datestr(gval(1),'HH:MM')];
   end;
   
elseif niceints(I)<=14,        % Ticks 1 day to 14 days apart
   labs=cellstr(datestr(gval+1/86400,'dd'));
   [Y,M,D,H,MI,S]=datevec(gval);
   kk=find(diff(M)~=0)+1;        % Pick out new months, label differently
   if any(kk),
     onelab=datestr(gval(kk),'dd-mmm-yyyy');
     onelab(onelab=='-')='/';
     labs(kk)=cellstr(onelab(:,[1:6]));
     labs{kk(1)}=onelab(1,[1:7 10 11]);
   else
     onelab=datestr(gval(1),'dd-mmm-yyyy');
     onelab(onelab=='-')='/';
     labs{1}=onelab([1:7 10 11]);;
   end;

elseif niceints(I)<365.25     % Ticks for various months
   labs=cellstr(datestr(gval+1/86400,'mmm'));
   [Y,M,D,H,MI,S]=datevec(gval);
   kk=find(diff(Y)~=0)+1;        % Pick out new years, label differently
   if any(kk),
     labs(kk)=cellstr(datestr(gval(kk),'yy'));
     labs(kk(1))=cellstr(datestr(gval(kk(1)),'yyyy'));
   else
     labs{1}=datestr(gval(1),'yyyy');
   end;
   
else      % Years
   labs=cellstr(datestr(gval,'yy'));
   labs{1}=datestr(gval(1),'yyyy');
end;

if nargout==0,
  set(gca,[ax 'tick'],gval,[ax 'ticklabel'],labs);
end;
