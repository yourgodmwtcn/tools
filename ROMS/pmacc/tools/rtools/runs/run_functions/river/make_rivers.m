function make_rivers(tspan,yearrun, S, gridfile, out_dir, riverFile, flowFile, initialize)
% make_rivers.m  12/6/2011  Parker MacCready
%
% edited by DAS, 2/9/2009
% sng updated to be run from run_maker.m
% sng updated to include a flag for river initialization 1/5/2012
% sng updated to include tspan and FlowFile such that you can input a flow
% files with multiple years in it and just extract the time span (tspan)
% that you desire July 2012
%
% tspan: time span that you want to extract flow data for
% yearrun: tag of year to run
% S: s-coordinate structure, same as in .in file
% gridfile: grid file for run with bathymetry
% out_dir: where to write rivers.nc to
% riverFile: it uses lat/lon to carve out rivers
%            in the current grid (if none, then need to set river variables
%            explicitly below) see Z_carve_river_channels and
%            Z_river_channels for information on what is in riverFile
% flowFile: a matlab file of river flow time series with the fields Qr_flow
%           Qr_year, Qr_yearday, T_riv, rname
% initialize is a flag to indicate whether you want the river flow to ramp
% up (1) or not (0) (only needed when starting a run from climatology)
%
%
% this creates a separate forcing file for time-dependent river flow.  It
% is slightly modified from code I got from John Warner at USGS
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Create a netcdf file that contains river forcing data for ROMS.
% Forcing data consists of:
% 'river_Xposition'  -   'river runoff  XI-positions at RHO-points'
% 'river_Eposition'  -   'river runoff ETA-positions at RHO-points'
% 'river_direction'  -   'river runoff direction'
% 'river_Vshape'     -   'river runoff mass transport vertical profile'
% 'river_transport'  -   'river runoff mass transport'
% 'river_flag'       -   'river runoff flag'
% 'river_temp'       -   'river runoff potential temperature'
% 'river_salt'       -   'river runoff salinity'
% 'river_mud_'       -   'river runoff suspended sediment concentration'
% 'river_sand_'      -   'river runoff suspended sediment concentration'
%
% In cppdefs.h you should have
% #define TS_PSOURCE
% #define UV_PSOURCE
% #undef  ANA_PSOURCE
% #undef  RIVER_SEDIMENT
%
%********************************************************
%  Begin user input section                             %
%********************************************************
tic

if nargin<7
    initialize = 1; %default is to initialize
end

%1) Create name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
riv_file_out=[out_dir 'rivers.nc'];

%1.5) specify the name of the river input data file, should have same number
%     of rivers as specified below or in riverFile in same order
load(flowFile);
% creates the variables (e.g. for 2004):
% note that rfile_in also requires these following variables plus it should
% have rname as a string of river names to keep things straight!
%   Qr_flow          nr x time                   double array of river flow (m^3/s)
%   Qr_year          1 x time                    double array (e.g. 2004)
%   Qr_yearday       1 x time                    double array (year day where e.g. 0 = Jan 1 2004)
%   T_riv            nr x time                   double array (degrees C)

% interpolate to fill in data at every day
Qr_time = Qr_yearday+datenum(Qr_year,1,1); %full Qr time in matlab date
Qr_yearday2 = Qr_time-datenum(str2double(yearrun),1,1); %year day in reference year
[~,num_rivers] = size(Qr_flow);
for i = 1:num_rivers
    mask = isnan(Qr_flow(:,i));
    disp(['Filling in ',num2str(sum(mask)),' NaNs in river: #' num2str(i) ' ' rname{i}]);
    Qr_flow(:,i) = interp1(Qr_yearday2(~mask),Qr_flow(~mask,i),Qr_yearday2);
    mask = isnan(T_riv(:,i));
    disp(['Filling in ',num2str(sum(mask)),' NaNs in Temp: #' num2str(i) ' ' rname{i}]);
    T_riv(:,i) = interp1(Qr_yearday2(~mask),T_riv(~mask,i),Qr_yearday2);
end


%2) Enter times of river forcings data, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
rtime = tspan(1):tspan(end);
river_time = (rtime-datenum(str2double(yearrun),1,1))*86400;
% would be zero at the start of the year
% the data is from USGS daily averages so having that be at 0:00 is maybe
% not perfect, but pretty good!
num_river_times=length(river_time);
% do not change this.

%5) Enter value of h, Lm, and Mm.
% Get some grid info, do not change this.
h=nc_varget(gridfile,'h');
[MP,LP]=size(h);

%********************************************************
% calc some grid stuff here - do not change this.
%********************************************************
Lm = LP-2;
Mm = MP-2;
N = S.N;
% **************************
L  = Lm+1;
M  = Mm+1;

% riverFile contains lat/lon positions in an array of structures named rivers
% where each structure has lat/lon, name, depth, and width
% *** note the last lat/lon of each river is where the point source should go.
load(riverFile);

% input riverFile to function Z_river_channels.m to alter bathymetry and
%   the mask to make sure river channels there
% NOTE that Z_river_channels does not actually alter the grid, to do
% that you need Z_carve_river_channels!!
rout = Z_river_channels(out_dir,rivers);

% now initialize river locations, directions, and flags
river_Xposition = rout(1).X(:);
river_Eposition = rout(1).Y(:);
river_direction = rout(1).D(:);
sign = rout(1).sign(:);
numcells = length(rout(1).X)*ones(length(rout(1).X),1);
river_ID = ones(length(rout(1).X),1);
numr = length(rivers);
if(numr > 1) %more than 1 river
    for j=2:numr
        river_Xposition = [river_Xposition;rout(j).X(:)];
        river_Eposition = [river_Eposition;rout(j).Y(:)];
        river_direction = [river_direction;rout(j).D(:)];
        sign = [sign;rout(j).sign(:)];
        numcells = [numcells; length(rout(j).X)*ones(length(rout(j).X),1)];
        river_ID = [river_ID; j*ones(length(rout(j).X),1)];
    end
end
Nsources = length(river_Xposition); % could be different than # of rivers
river_flag = 3*ones(Nsources,1); % 3 specifies temp and sal at each point source

%8) Initialize river shape.
for i=1:Nsources
    river_Vshape(:,i) = flipud(linspace(2/N,0,N)');
    %linear decay form surface value to 0, fixed sng 7/2011
end

%9) Initialize river flow.
river_transport = nan(num_river_times,Nsources);
for i=1:Nsources
    fac = ones(size(river_time));
    if initialize==1
        fac(river_time<43200) = 1.0+tanh((river_time(time)-43200.0)/43200.0);
        disp('ramping up river flow');
    end
    Qi = interp1(Qr_time,Qr_flow(:,river_ID(i)),rtime); %interpolate to the desired time span
    Q = Qi./numcells(i); %divide this river Q by number of cells
    river_transport(:,i)=sign(i).*fac.*Q; %make negative for S or W
end

%10) Time series of river temp and salt.
river_temp = nan(num_river_times,N,Nsources);
river_salt = zeros(num_river_times,N,Nsources);
for k=1:N
    for i=1:Nsources
        river_temp(:,k,i) = interp1(Qr_time,T_riv(:,river_ID(i)),rtime);
    end
end

%
NAT=2;  %assume temp + salt are active
NT = NAT;        % total number of tracers.

%********************************************************
%  END of USER INPUT
%********************************************************

my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
disp(['*** creating ',riv_file_out,' ***']);
nc_create_empty ( riv_file_out, my_mode );
nc_padheader (riv_file_out , 30000 );

% Global attributes:
disp(' ## Defining Global Attributes...')
nc_attput(riv_file_out,nc_global,'history', ...
    ['Created by "' mfilename '" on ' datestr(now)]);
nc_attput(riv_file_out,nc_global,'type', ...
    'Initialization file from create_roms_init.m');

% Dimensions:
disp(' ## Defining Dimensions...')

nc_add_dimension(riv_file_out, 'xi_rho', LP);
nc_add_dimension(riv_file_out, 'eta_rho', MP);
nc_add_dimension(riv_file_out, 'xi_psi', L);
nc_add_dimension(riv_file_out, 'eta_psi', M);
nc_add_dimension(riv_file_out, 'xi_u', L);
nc_add_dimension(riv_file_out, 'eta_u', MP);
nc_add_dimension(riv_file_out, 'xi_v', LP);
nc_add_dimension(riv_file_out, 'eta_v', M);

nc_add_dimension(riv_file_out, 's_rho', N);
nc_add_dimension(riv_file_out, 'tracer', NT);
nc_add_dimension(riv_file_out, 's_w', N+1);

nc_add_dimension(riv_file_out, 'one', 1);
nc_add_dimension(riv_file_out, 'two', 2);
nc_add_dimension(riv_file_out, 'river', Nsources);
nc_add_dimension(riv_file_out, 'river_time', num_river_times);

% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

varstruct.Name = 'theta_b';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate bottom control parameter'];
units = '1';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'theta_s';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate surface control parameter'];
units = '1';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Tcline';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate surface/bottom layer width'];
units = 'meter';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'hc';
varstruct.Dimension = {'one'};
long_name = ['S-coordinate parameter, critical depth'];
units = 'meter';
varstruct.Attribute = struct('Name', ...
    {'long_name','units'},'Value',{long_name,units});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Cs_r';
varstruct.Dimension = {'s_rho'};
long_name = ['S-coordinate stretching curves at RHO-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'Cs_r, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'Cs_w';
varstruct.Dimension = {'s_w'};
long_name = ['S-coordinate stretching curves at W-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'Cs_w, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'sc_r';
varstruct.Dimension = {'s_rho'};
long_name = ['S-coordinate at RHO-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'sc_r, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'sc_w';
varstruct.Dimension = {'s_w'};
long_name = ['S-coordinate at W-points'];
units = '1';
valid_min = -1;
valid_max = 0;
field = 'sc_w, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','valid_min','valid_max','field'},...
    'Value',{long_name,units,valid_min,valid_max,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river';
varstruct.Dimension = {'river'};
long_name = ['river_runoff identification number'];
units = 'nondimensional';
field = 'num_rivers, scalar';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_time';
varstruct.Dimension = {'river_time'};
long_name = ['river_time'];
units = 'seconds';
field = 'river_time, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Xposition';
varstruct.Dimension = {'river'};
long_name = ['river runoff  XI-positions at RHO-points'];
units = 'scalar';
field = 'river runoff XI position, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Eposition';
varstruct.Dimension = {'river'};
long_name = ['river runoff  ETA-positions at RHO-points'];
units = 'scalar';
field = 'river runoff ETA position, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_direction';
varstruct.Dimension = {'river'};
long_name = ['river runoff direction, XI=0, ETA>0'];
units = 'scalar';
field = 'river runoff direction, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_Vshape';
varstruct.Dimension = {'s_rho','river'};
long_name = ['river runoff mass transport vertical profile'];
units = 'scalar';
field = 'river runoff vertical profile, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_transport';
varstruct.Dimension = {'river_time','river'};
long_name = ['river runoff mass transport'];
units = 'meter^3/s';
field = 'river runoff mass transport, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_flag';
varstruct.Dimension = {'river'};
long_name = ['river flag, 1=temp, 2=salt, 3=temp+salt, 4=temp+salt+sed, 5=temp+salt+sed+bio'];
units = 'nondimensional';
field = 'river flag, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_temp';
varstruct.Dimension = {'river_time','s_rho','river'};
long_name = ['river runoff potential temperature'];
units = 'Celsius';
field = 'river temp, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

varstruct.Name = 'river_salt';
varstruct.Dimension = {'river_time','s_rho','river'};
long_name = ['river runoff salinity'];
units = 'psu';
field = 'river salinity, scalar, series';
varstruct.Attribute = struct('Name', ...
    {'long_name','units','field'},...
    'Value',{long_name,units,field});
nc_addvar(riv_file_out, varstruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')
nc_varput(riv_file_out, 'theta_s', S.theta_s);
nc_varput(riv_file_out, 'theta_b', S.theta_b);
nc_varput(riv_file_out, 'Tcline', S.tcline);
nc_varput(riv_file_out, 'Cs_r', S.Cs_r);
nc_varput(riv_file_out, 'Cs_w', S.Cs_w);
nc_varput(riv_file_out, 'sc_w', S.s_w);
nc_varput(riv_file_out, 'sc_r', S.s_rho);
nc_varput(riv_file_out, 'hc', S.hc);
nc_varput(riv_file_out, 'river', river_ID);

nc_varput(riv_file_out, 'river_time', river_time);
nc_varput(riv_file_out, 'river_Xposition', river_Xposition);
nc_varput(riv_file_out, 'river_Eposition', river_Eposition);
nc_varput(riv_file_out, 'river_direction', river_direction);
nc_varput(riv_file_out, 'river_Vshape', river_Vshape);
nc_varput(riv_file_out, 'river_transport', river_transport);
nc_varput(riv_file_out, 'river_flag', river_flag);
nc_varput(riv_file_out, 'river_temp', river_temp);
nc_varput(riv_file_out, 'river_salt', river_salt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('DONE')
trivend = toc;
disp(['rivers file complete in ',num2str(round(trivend)),' sec']);



