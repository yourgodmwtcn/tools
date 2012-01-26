function mm = mm_setup(varargin)
% MM_SETUP Set up a new movieman instance.
%   MM = MM_SETUP('param1','value1','param2','value2', etc. ... )
%   Creates a movieman structure MM with the specified parameters.
%   All of the parameters are optional. The possible parameteris,
%   along with their default values are listed below:
% 
%   outputFile: 'mm_output.mpg' Filename for the output file.
%
%     frameDir: './mm_frames/'  Temporary directory for storing frames.
%
%      bitRate: 25000000        Specifies the bits/sec (determines quality).
%
%    frameRate: 25              The frame rate (frames/sec).
%                               Only 25 and 30 are allowed for mpeg-2 video.
%
%    pixelSize: [800 300]       Dimensions of the movie in pixels [x y].
%
%   ffmpegArgs: ''              Any extra arguments you might want to pass
%                               ffmpeg. (Be careful with this, advanced
%                               users only.)
%
%   For more info see help for mm_addFrame and mm_render.
%
%   Copyright (c) 2009 Ryan Abernathey
%   http://web.mit.edu/rpa/www/movieman
%   See LICENSE file for license information. 
%

p = inputParser;
p.addParamValue('frameDir','./mm_frames/',@check_frameDir);
p.addParamValue('pixelSize',[800 300],@check_pixelSize);
p.addParamValue('outputFile','mm_output.mpg');
p.addParamValue('frameRate',25,@isnumeric);
p.addParamValue('bitRate',25000000,@isnumeric);
p.addParamValue('ffmpegArgs','');

p.parse(varargin{:});
mm = struct(p.Results);

% check if we need to make the frame dir (using defaults)
if any(strcmpi('frameDir',p.UsingDefaults))
    check_frameDir(mm.frameDir);
end
% check the output directory
files = dir([mm.frameDir '/mm_frame_*']);
if (length(files)>0)
    ('Found pre-existing frames in the frameDir. This could confuse MovieMan.');
    reply = input('Do you want to delete the frames? Y/N [Y]: ', 's');
    if isempty(reply)
        reply = 'Y';
    end
    if strcmpi(reply,'Y')
        delete([mm.frameDir '/mm_frame_*']);
    end        
end

end

function ok = check_frameDir(d)
    if (~ischar(d))
        error('frameDir argument must be a char');
    end
    if (isdir(d))
        disp(['frameDir "' d '" already exists (that''s okay!)']);
        ok=true;
    else
        disp(['frameDir "' d '" doesn''t exist, creating it now']);
        ok=mkdir(d);
    end    
end

function ok = check_pixelSize(p)
    if (length(p) ~= 2)
        ok=false;
    else
        ok=all(mod(p,1)==0);
    end
end
