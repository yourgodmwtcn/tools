%       [] = roms_browse(fname, varname, tindices, axis, index)

function [] = roms_browse(fname, varname, tindices, axis, index)

ckey = ''
load trees
 figure;
 imshow(X,map);
 tmp = get(gcf,'currentkey');
 while strcmp(ckey,tmp)
 disp('Waiting for keypress!');
 pause(0.1);
 ckey = get(gcf,'currentkey');
 end
 
 ckey
    