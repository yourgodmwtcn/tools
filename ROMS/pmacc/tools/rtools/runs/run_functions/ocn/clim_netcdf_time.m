function [] = clim_netcdf_time(outfile,V,ntimes)
% 12/6/2011  Parker MacCready
% adds the time dimension, if appropriate

% define rest of variables and attributes
if V.do_addtimevar
    nc_add_dimension(outfile, V.nctimename, ntimes);
    varstruct.Name = V.nctimename;
    varstruct.Dimension = {V.nctimename};
    long_name = [V.nclongname,' time'];
    units = 'seconds';
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(outfile,varstruct);
end

