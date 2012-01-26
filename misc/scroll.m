function scroll(old_val)
%  scroll(old_val)
%
% SCROLL - Scroll subplots vertically
%   Used by scrollsubplot.
%
%   See also SCROLLSUBPLOT

% Copyright Bjorn Gustavsson 20050526


fig = gcf;
% Get all figure's children
fig_chld = get(gcf,'children');
% And their type
for i =1:length(fig_chld),    
  chld_type{i} = get(fig_chld(i),'Type');
end

% Look among the 'callback' ones
clbk_indx = strmatch('uicontrol',chld_type);

clbk = [];
for i = clbk_indx
  clbk{end+1} = get(fig_chld(i),'callback');
end

if ~isempty(clbk)
  scroll_clbk_indx = strmatch('scroll',clbk);
end


if isempty(scroll_clbk_indx)
  return
end
clbk_ui = fig_chld(clbk_indx(scroll_clbk_indx));

a_indx = strmatch('axes',chld_type);

for i = a_indx,
  try 
    a_pos(i,:) = cell2mat([get(fig_chld(i),'position')]);
  catch
    a_pos(i,:) = get(fig_chld(i),'position');
  end
  
end

pos_y_range = [min(.07,min(a_pos(a_indx,2))) max(a_pos(a_indx,2) + a_pos(a_indx,4) )+.07-.9];

val = get(clbk_ui,'value');
[ val  old_val  diff(pos_y_range);];
step = ( old_val - val) * diff(pos_y_range);


for i = 1:length(a_indx),
  
  set(fig_chld(a_indx(i)),'position',get(fig_chld(a_indx(i)),'position') + [0 step 0 0]);
  
end

old_val = val;

set(fig_chld(clbk_indx(scroll_clbk_indx)),'callback',['scroll(',num2str(val),')']);
