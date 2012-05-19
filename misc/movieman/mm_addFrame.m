function ok = ff_addFrame(mm,f,n)
% MM_ADDFRAME Add a frame to your movieman movie.
%   FF_ADDFRAME(MM,FIGNUM)
%    MM is a movieman struct (created using MM_SETUP).
%    FIGNUM is a figure handle.
%
%   FF_ADDFRAME(MM,FIGNUM,N)
%    The optional argument N prints the frame N times.
%
%   For more info see help for mm_setup and mm_render.
%
%   Copyright (c) 2009 Ryan Abernathey
%   http://web.mit.edu/rpa/www/movieman
%   See LICENSE file for license information. 
%


if isempty(mm.frameDir)
    error('MovieMan struct mm is incomplete!');
end

if nargin < 3, n = 1; end

% set up the resolution
dpi = 96;
dpispec = sprintf('-r%02d',dpi);
inchsize = mm.pixelSize / dpi;


set(f,'Units','inches');
set(f,'PaperPositionMode','manual','PaperSize',inchsize,'PaperPosition',[0 0 inchsize]);
set(f,'InvertHardcopy','off')
set(f,'PaperPositionMode','manual','PaperSize',inchsize,...
    'PaperPosition',[0 0 inchsize]);

% figure out the frame number from the directory
frames = dir([mm.frameDir '/mm_frame_*']);
if (length(frames)==0)
    N = 0;
else
    lastframe = frames(end).name;
    N = str2num(lastframe(10:15));
end

% output it
for i=1:n    
    filename = sprintf('mm_frame_%06d',N+i);
    % for PPM
    %print(dpispec,'-zbuffer','-dppmraw', [mm.frameDir '/' filename '.ppm']);
    % for PNG
    print(dpispec,'-zbuffer','-dpng', [mm.frameDir '/' filename '.png']);
    %export_fig(dpispec,[mm.frameDir '/' filename '.png']);
    fprintf('Just added frame %6d\n',N+i);
end

ok = true;

end