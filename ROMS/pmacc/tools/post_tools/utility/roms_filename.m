function fname = roms_filename(basename, n);

% fname = roms_filename(basename, n);

nstr = ['0000' num2str(n)];
nstr = nstr(end-3:end);
fname = [basename nstr '.nc'];
