function theAxis = scrollsubplot(nrows, ncols, thisPlot)
%SCROLLSUBPLOT Create axes in tiled positions.
%   SCOLLSUBPLOT(m,n,p), breaks the Figure window into
%   an m-by-n matrix of small axes, selects the p-th axes for 
%   for the current plot, and returns the axis handle.  The axes 
%   are counted along the top row of the Figure window, then the
%   second row, etc.  For example,
% 
%       SCROLLSUBPLOT(3,1,-1), PLOT(income)
%       SCROLLSUBPLOT(3,1,1), PLOT(declared_income)
%       SCROLLSUBPLOT(3,1,2), PLOT(tax)
%       SCROLLSUBPLOT(3,1,3), PLOT(net_income)
%       SCROLLSUBPLOT(3,1,4), PLOT(the_little_extra)
%   plots declared_income on the top third of the window, tax in
%   the middle, and the net_income in the bottom third. Above the
%   top of the figure income is ploted and below the lower edge
%   the_little_extra is to be found. To navigate there is a slider
%   along the right figure edge.
% 
%   The function works well for regular grids where m,n is constant
%   for all p. When m,n varies there is no guaranti that the steps
%   of the slider is nicely adjusted to the sizes of the
%   subplots.
%
%   Differences with SUBPLOT: SCROLLSUBPLOT requires 3 input
%   arguments, no compatibility with subplot(323), no handle as
%   input. Further  PERC_OFFSET_L is decreased from 2*0.09 to 0.07
%   and PERC_OFFSET_R is decreased from 2*0.045 to 0.085. This
%   leaves less space for titles and labels, but give a plaid grid
%   of subplots even outside the visible figure area.
%   
%   Bug/feature when the slider is shifted from its initial
%   position and then extra subplots is added, they get
%   mis-positioned.
%   
%   See also SCROLL, SUBPLOT,

%   Copyright Bjorn Gustavsson 20050526, Modification/extension of
%   Mathworks subplot. GPL license apply.


persistent maxrownr minrownr

%%% This is a matlab bug(?) that blanks axes that have been rescaled to
%%% accomodate a colorbar. I've not tried to fix it. BG
% we will kill all overlapping siblings if we encounter the mnp
% specifier, else we won't bother to check:
narg = nargin;
kill_siblings = 0;
create_axis = 1;
delay_destroy = 0;
tol = sqrt(eps);
if narg ~= 3 % not compatible with 3.5, i.e. subplot ==
             % subplot(111) errors out
  error('Wrong number of arguments')
end

%check for encoded format
handle = '';
position = '';

kill_siblings = 1;

% if we recovered an identifier earlier, use it:
if(~isempty(handle))
  
  set(get(0,'CurrentFigure'),'CurrentAxes',handle);
  
elseif(isempty(position))
  % if we haven't recovered position yet, generate it from mnp info:
  if (min(thisPlot) < 1)&0
    error('Illegal plot number.')
  else
    % This is the percent offset from the subplot grid of the plotbox.
    PERC_OFFSET_L = 0.07;
    PERC_OFFSET_R = 0.085;
    PERC_OFFSET_B = PERC_OFFSET_L;
    PERC_OFFSET_T = PERC_OFFSET_R;
    if nrows > 2
      PERC_OFFSET_T = 0.9*PERC_OFFSET_T;
      PERC_OFFSET_B = 0.9*PERC_OFFSET_B;
    end
    if ncols > 2
      PERC_OFFSET_L = 0.9*PERC_OFFSET_L;
      PERC_OFFSET_R = 0.9*PERC_OFFSET_R;
    end

    % Subplots version:
    % row = (nrows-1) -fix((thisPlot-1)/ncols)
    % col = rem (thisPlot-1, ncols);
    % Slightly midified to alow for having ngative thisPlot
    row = (nrows-1) -floor((thisPlot-1)/ncols);
    col = mod (thisPlot-1, ncols);
    
    % From here on to line 190 essentially identical to SUBPLOT (==untouched)
    
    % For this to work the default axes position must be in normalized coordinates
    if ~strcmp(get(gcf,'defaultaxesunits'),'normalized')
      warning('DefaultAxesUnits not normalized.')
      tmp = axes;
      set(axes,'units','normalized')
      def_pos = get(tmp,'position');
      delete(tmp)
    else
      def_pos = get(gcf,'DefaultAxesPosition')+[-.05 -.05 +.1 +.05];
    end
    col_offset = def_pos(3)*(PERC_OFFSET_L+PERC_OFFSET_R)/ ...
	(ncols-PERC_OFFSET_L-PERC_OFFSET_R);
    row_offset = def_pos(4)*(PERC_OFFSET_B+PERC_OFFSET_T)/ ...
	(nrows-PERC_OFFSET_B-PERC_OFFSET_T);
    totalwidth = def_pos(3) + col_offset;
    totalheight = def_pos(4) + row_offset;
    width = totalwidth/ncols*(max(col)-min(col)+1)-col_offset;
    height = totalheight/nrows*(max(row)-min(row)+1)-row_offset;
    position = [def_pos(1)+min(col)*totalwidth/ncols ...
		def_pos(2)+min(row)*totalheight/nrows ...
		width height];
    if width <= 0.5*totalwidth/ncols
      position(1) = def_pos(1)+min(col)*(def_pos(3)/ncols);
      position(3) = 0.7*(def_pos(3)/ncols)*(max(col)-min(col)+1);
    end
    if height <= 0.5*totalheight/nrows
      position(2) = def_pos(2)+min(row)*(def_pos(4)/nrows);
      position(4) = 0.7*(def_pos(4)/nrows)*(max(row)-min(row)+1);
    end
  end
end

% kill overlapping siblings if mnp specifier was used:
nextstate = get(gcf,'nextplot');
if strncmp(nextstate,'replace',7), nextstate = 'add'; end
if(kill_siblings)
  if delay_destroy
    if nargout
      error('Function called with too many output arguments')
    else
      set(gcf,'NextPlot','replace'); return,
    end
  end
  sibs = get(gcf, 'Children');
  got_one = 0;
  for i = 1:length(sibs)
    if(strcmp(get(sibs(i),'Type'),'axes'))
      units = get(sibs(i),'Units');
      set(sibs(i),'Units','normalized')
      sibpos = get(sibs(i),'Position');
      set(sibs(i),'Units',units);
      intersect = 1;
      if(     (position(1) >= sibpos(1) + sibpos(3)-tol) | ...
	      (sibpos(1) >= position(1) + position(3)-tol) | ...
	      (position(2) >= sibpos(2) + sibpos(4)-tol) | ...
	      (sibpos(2) >= position(2) + position(4)-tol))
	intersect = 0;
      end
      if intersect
	if got_one | any(abs(sibpos - position) > tol)
	  delete(sibs(i));
	else
	  got_one = 1;
	  set(gcf,'CurrentAxes',sibs(i));
	  if strcmp(nextstate,'new')
	    create_axis = 1;
	  else
	    create_axis = 0;
	  end
	end
      end
    end
  end
  set(gcf,'NextPlot',nextstate);
end

% create the axis:
if create_axis
  if strcmp(nextstate,'new'), figure, end
  ax = axes('units','normal','Position', position);
  set(ax,'units',get(gcf,'defaultaxesunits'))
else 
  ax = gca; 
end


% return identifier, if requested:
if(nargout > 0)
  theAxis = ax;
end



%%% From here on out set up scrollbar if needed
scroll_clbk_indx = [];

% Get all figure's choldren
fig_chld = get(gcf,'children');
% And their type
for i =1:length(fig_chld),    
  chld_type{i} = get(fig_chld(i),'Type');
end
%chld_type = chld_type
% Look among the 'callback' ones
clbk_indx = strmatch('uicontrol',chld_type);

ax_indx = strmatch('axes',chld_type);
if length(ax_indx)==1
  maxrownr = -inf;
  minrownr = inf;
end

clbk = [];
for i = clbk_indx
  clbk{end+1} = get(fig_chld(i),'callback');
end

if ~isempty(clbk)
  scroll_clbk_indx = strmatch('scroll',clbk);
end

if isempty(scroll_clbk_indx)
  uicontrol('Units','normalized',...
            'Style','Slider',...
            'Position',[.98,0,.02,1],...
            'Min',0,...
            'Max',1,...
            'Value',1,...
            'visible','off',...
            'Callback','scroll(1)');
end
if ( nrows*ncols < thisPlot | thisPlot < 1 )
  set(fig_chld(clbk_indx(scroll_clbk_indx)),'visible','on')
end
scroll(1)

maxrownr = max(maxrownr,-row);
lminrnr = minrownr;
minrownr = min(minrownr,-row);
set(fig_chld(clbk_indx(scroll_clbk_indx)),...
    'sliderstep',[1/nrows 1]/(1/((nrows)/max(1,1+maxrownr-minrownr-nrows))))
set(fig_chld(clbk_indx(scroll_clbk_indx)),...
    'value',1)
