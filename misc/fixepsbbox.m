function fixepsbbox(filename)
% function fixepsbbox(filename)
% 
% matlab seems to compute a bounding box on eps files which is too
% large in the x-direction
% 
% this script fixes the bounding box
% it is 99% stolen from fixeps.m, located here:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4818&objectType=file
% the only change is that this changes the bbox numbers
% seriously, the only change is the addition of lines 22,23,33,34
% 
% boundingbox has form of:
% %%BoundingBox:    x1   y1   x2   y2
% where (x1,y1) is lower-left and (x2,y2) is upper-right
% 
% matlab computes x1 too small and x2 too large
% changes lines to:
% %%BoundingBox:    x1+dx1 y1 x2+dx2 y2

% default amount to change bbox - found this fixed my plots just fine
dx1 = 10;	% amount to move x1
dx2 = -25;	% amount to move x2
dy1 = -10;
dy2 = -10;

fid = fopen(filename,'r+');
k=0;
while k <2                                                  % 2 locations to replace.
    tline = fgetl(fid);                                     % get one line text
    stridx=strfind(tline,'Box:');
	if isempty(stridx)==0
        len=length(tline);                                  % the original line length
		bb=sscanf(tline(stridx+4:end),'%i');                % read the numbers
		bb(1) = bb(1) + dx1;% change x1
        bb(2) = bb(2) + dy1;
		bb(3) = bb(3) + dx2;								% change x2
        bb(4) = bb(4) + dy2;
		bbstr=sprintf('%g %g %g %g',bb);                    % write bb numbers to string
        tline=tline(1:stridx+3);                            % keep the "%%(page)boundingbox" string (with starting '%%')
		spaces(1:len-length(tline)-length(bbstr)-1)=' ';    % add trailing spaces as to overwrite old line completely
		tline=[tline ' ' bbstr spaces];                     % concate numbers and blank spaces to "%%(page)boundingbox"

        fseek(fid,-len-2,'cof');                            % before using fprintf search to correct position
		count = fprintf(fid,'%s',tline);
        fseek(fid,2,'cof');                                 % seek to beginning of line (for windows text file) on
                                                            % for linux: change '2' to '1' I think
        k=k+1;
	end
end
fclose(fid);
