function [] = make_atm(gridfile,year,dth,indir0,outdir)
% 12/18/2012  Parker MacCready
%
% creates ROMS surface forcing fields from our heavily-preprocessed MM5
% and WRF files

year0 = str2num(year);

% make the output directory if needed
if ~exist(outdir,'dir'); mkdir(outdir); end;

% here is vector of times we want,
% including some padding into the next year
tdvec = [datenum(year0,1,1,0,0,0):dth/24: ...
    datenum(year0+1,1,1,0,0,0)+(dth/24)]';

% list of input variables
varlist = {'psfc';'t2';'qair';'rain';'swdown';'lwdown';'u10r';'v10r'};

% create index vectors to help find things (note that the processed atm
% files are stored in month files)
% it doesn't matter which variable we choose
var = varlist{1};
%
td = []; % initialize datenum vector
mo = []; % initialize a month vector
yr = []; % initialize a year vector
ind = []; % index within the month
data = []; % initialize a data vector
flag = []; % intialize flag for interpolated fields
%
for year = year0:year0+1
    indir = [indir0,num2str(year),'/'];
    for month = 1:12
        load([indir,var,'_',num2str(month),'.mat'],'TD');
        nt = length(TD);
        td = [td; TD];
        mo = [mo; month*ones(nt,1)];
        yr = [yr; year*ones(nt,1)];
        ind = [ind; [1:nt]'];
    end % end of month loop
end % end of year loop

% Read in the grid
lon_rho = nc_varget(gridfile, 'lon_rho');
lat_rho = nc_varget(gridfile, 'lat_rho');
[M,L] = size(lon_rho);

for vv = 1:length(varlist);
    tic
    var = varlist{vv};
    
    disp(['* Writing forcing file for ',var]);
    % *** create and populate the NetCDF file for ROMS ***
    
    % vartime units = seconds from the start of the year
    vartime = round(86400*(tdvec - datenum(year0,1,1,0,0,0)));
    
    scalefactor = 1; % multiply by this
    scalefactor2 = 0; % and then add this
    switch var
        case 'psfc'
            ncvarname = 'Pair';
            nclongname = 'surface air pressure';
            ncunits = 'millibar';
            nctimename = 'pair_time';
            scalefactor = 1/100; % convert Pa to mbar
        case 'rain'
            ncvarname = 'rain';
            nclongname = 'rain fall rate';
            ncunits = 'kilograms meter-2 second-2';
            nctimename = 'rain_time';
            scalefactor = 10; % convert cm s-1 to kg m-2 s-1
        case 'swdown'
            ncvarname = 'swrad';
            nclongname = 'solar shortwave radiation flux';
            ncunits = 'watts meter-2';
            nctimename = 'srf_time';
            scalefactor = 1 - 0.1446; % account for reflection
        case 'lwdown'
            ncvarname = 'lwrad_down';
            nclongname = 'downwelling longwave radiation flux';
            ncunits = 'watts meter-2';
            nctimename = 'lrf_time';
        case 't2'
            ncvarname = 'Tair';
            nclongname = 'surface air temperature';
            ncunits = 'Celsius';
            nctimename = 'tair_time';
            scalefactor2 = -273.15; % convert K to C
        case 'qair'
            ncvarname = 'Qair';
            nclongname = 'surface air relative humidity';
            ncunits = 'percentage';
            nctimename = 'qair_time';
        case 'u10r'
            ncvarname = 'Uwind';
            nclongname = 'surface u-wind component';
            ncunits = 'meter second-1';
            nctimename = 'wind_time';
        case 'v10r'
            ncvarname = 'Vwind';
            nclongname = 'surface v-wind component';
            ncunits = 'meter second-1';
            nctimename = 'wind_time';
    end
    
    % create the NetCDF file
    frcname = [outdir, ncvarname, '.nc'];
    
    disp('  ++ Creating the NetCDF file')
    my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
    nc_create_empty( frcname, my_mode );
    % global attributes
    nc_padheader (frcname, 20000 );
    nc_attput(frcname, nc_global, 'type','ROMS Surface Forcing File');
    
    % define dimensions
    nc_add_dimension(frcname, 'xi_rho', L);
    nc_add_dimension(frcname, 'eta_rho', M);
    nc_add_dimension(frcname, nctimename, length(vartime));
    
    % define variables and attributes
    %   first time
    varstruct.Name = nctimename;
    varstruct.Dimension = {nctimename};
    long_name = [nclongname,' time'];
    units = 'seconds';
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    %   then the field
    varstruct.Name = ncvarname;
    varstruct.Dimension = {nctimename, 'eta_rho', 'xi_rho'};
    long_name = nclongname;
    units = ncunits;
    varstruct.Attribute = struct('Name', ...
        {'long_name','units'},'Value',{long_name,units});
    nc_addvar(frcname, varstruct);
    
    % write the variables
    disp('  ++ Writing to the file')
    nc_varput(frcname, nctimename, vartime);
    
    %BEGIN TIME LOOP
    infile_prev = ['blank'];
    for iii = 1:length(vartime)
        
        % specify the input file
        this_td = tdvec(iii);
        
        % this call to dsearchn gets the closest time
        itd = dsearchn(td,this_td);
                
        indir = [indir0,num2str(yr(itd)),'/'];
        
        infile = [indir,var,'_',num2str(mo(itd)),'.mat'];
        
        % only reload if the name of indir changes
        if ~strcmp(infile,infile_prev)
            disp(['   loading (',num2str(yr(itd)),') ', ...
                var,'_',num2str(mo(itd)),'.mat'])
            load(infile);
            infile_prev = infile;
        end
        
        % VAR exists because it is in "infile"
        AA = squeeze(VAR(ind(itd),:,:));
        
        AAA = interp2(LON,LAT,AA,lon_rho,lat_rho);
        
        nc_varput(frcname, ncvarname, scalefactor*AAA + scalefactor2, ...
            [iii-1 0 0], [1 M L]);
        
    end % end of time loop
    
    dt_toc = toc;
    disp(['  took ',num2str(round(dt_toc)),' sec']);
    
end % end of variable loop

