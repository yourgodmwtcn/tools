function RD = roms_createRunDef(shortname, dirname, hisbasename, diabasename)

% RD = roms_createRunDef(dirname);
% RD = roms_createRunDef(shortname, dirname);
% RD = roms_createRunDef(shortname, dirname, hisbasename, diabasename);
%
% makes a structure containing the run shortname, dirname, seriesDefs
% for the _his and _dia files (if any _his or _dia files were found),
% and model grid
%
% neil banas feb 2009
% sng included option for nargin==3 april 2011

if nargin==1
	dirname = shortname;
	shortname = '';
	hisbasename = 'ocean_his_';
	diabasename = 'ocean_dia_';
end
if nargin==2
	hisbasename = 'ocean_his_';
	diabasename = 'ocean_dia_';
end
if nargin==3
    diabasename = 'ocean_dia_';
end

RD.shortname = shortname;
if dirname(end) ~= '/', dirname = [dirname '/']; end
RD.dirname = dirname;

his = roms_createSeriesDef(dirname, hisbasename);
if ~isempty(his.ncn)
	RD.his = his;
	RD.grid = roms_loadGrid(RD.his);
end

dia = roms_createSeriesDef(dirname, diabasename);
if ~isempty(dia.ncn)
	RD.dia = dia;
	if length(dia.nctime) > 1
		RD.dia.dt = RD.dia.nctime(2) - RD.dia.nctime(1); % time between saves: needed to convert the N_Flux_ terms into actual rates
	end
	if ~isfield(RD,'grid')  % if the grid wasn't loaded from the _his series
		RD.grid = roms_loadGrid(RD.dia);
	end
end