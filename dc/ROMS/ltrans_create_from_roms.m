% creates LTRANS initial location csv file based on ROMS float output
%   [] = ltrans_create_from_roms(fname,flt,rgrid)
% fname = output file name; flt = ROMS float file, rgrid

function [] = ltrans_create_from_roms(fname,flt,rgrid)

    fid = fopen(fname,'w');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    roms = floats('roms',flt,rgrid);
    
    fprintf(fid,'%.2f,%.2f,%.2f,%d \n',roms.init');
    fprintf('Added %d floats to %s \n',size(roms.x,2),fname);

    fclose(fid);