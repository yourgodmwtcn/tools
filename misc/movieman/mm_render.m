function mm_render(mm)
% MM_RENDER Render a movieman movie using ffmpeg.
%   FF_RENDER(MM)
%    Take a movieman struct MM and render a movie using the params
%    specified when MM was created (using mm_setup).
%
%   This function will not work if you don't have ffmpeg in your path.
%   (http://www.ffmpeg.org/)
%
%   For more info see help for mm_setup and mm_render.
%
%   Copyright (c) 2009 Ryan Abernathey
%   http://web.mit.edu/rpa/www/movieman
%   See LICENSE file for license information. 
%

% We don't specify the movie size because ffmpeg will automatically detect it
% based on the images size of the frames.
cmd = sprintf('ffmpeg -g %d -b %d -f image2 -i %s %s %s',...
    mm.frameRate,mm.bitRate,[mm.frameDir '/mm_frame_%06d.png'],...
    mm.ffmpegArgs,mm.outputFile);
%     mm.frameRate,mm.bitRate,[mm.frameDir '/mm_frame_%06d.ppm'],...

disp(cmd);

disp('Looking for ffmpeg...');
foundit = system('which ffmpeg');
if (foundit ~= 0)
    warning('The "ffmpeg" command was not found in your path.');
    disp('For information on installing ffmpeg, visit http://www.ffmpeg.org/.');
    error('Unable to render the movie without ffmpeg. Giving up.')
end

system(cmd);

end