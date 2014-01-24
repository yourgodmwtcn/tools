function EX2D = obs_mooringExtract(filename, vars, timerange)
%--------------------------------------------------------------------------
%       EX2D = obs_mooringExtract(filename, vars, timeRange)
%                                         
% filename must be string with only one filename
% looks inside a netcdf file for one or more observational variables and
% returns them in a structure _EX2D_, along with coordinate variables.
%
% vars can be a string containing one variable name, or a cell array of strings,
% or 'all'.
%
% timeRange  is a 2 element vector of matlab datenums or 'all'.
%

%
% written by C. Bassin ,D. Sutherland and N. Banas, UW, Jan 2010
%
%%

% if (nargin~=7 ||nargin~=6|| nargin~=4)
%     error ('Not correct number of arguments')
% end

% check filename exists:
if exist(filename, 'file')~=2
    warning(['Filename ' filename ' not found, returning empty structure']);
    EX2D=empty_structure_return(vars);
    return

elseif exist(filename, 'file')==2

    %check variables exist in file
    %first check that time, lat and lon and pressure or depth exist which are necessary for
    %program to run prperly
    if nc_isvar(filename,'time')==0;
        EX2D=empty_structure_return(vars);
        return
    end
    if nc_isvar(filename,'longitude')==0;
        EX2D=empty_structure_return(vars);
        return
    end
    if nc_isvar(filename,'latitude')==0;
        EX2D=empty_structure_return(vars);
        return
    end
    if nc_isvar(filename,'pressure')==0
        if nc_isvar(filename,'depth')==0
            EX2D=empty_structure_return(vars);
            return
        end
    end
    % how many variables, check if it is a variable name or cell array
       if strcmp(vars,'all')
           Info_file=nc_info(filename);
           for v=1:length(Info_file.Dataset)
               vars_in_file{v}=Info_file.Dataset(v).Name;
           end
           vars=vars_in_file;
       end
    
    
    if ischar(vars)==1;
        var_exist=nc_isvar(filename,vars);

    elseif iscell(vars)==1
        %check whether individual variables exist

        var_exist=nan(length(vars),1);
        for L=1:length(vars);
            var_exist(L)=nc_isvar(filename,vars{L});
        end
    end

    %no variables EX2Dist in files
    if sum(var_exist)==0
        % warning(['no requested variables exist in the files' filename ', returning nan structure']);
        EX2D=empty_structure_return(vars);
        return
    end




    %{
                at least one variable exists, so get lat, lon, associate
                within time frame, then get variables that exist
    %}

    longitude=nc_varget(filename,'longitude');
    latitude=nc_varget(filename,'latitude');
   

 
    %get time array and check that time exist

    T=nc_varget(filename,'time');  % use mean time to cheeck columns in case a profile goes into the next day
    if ischar(timerange)==1
        if strcmpi(timerange,'all')==1
            TimeInFile=1:length(T);
        end
    else

        TimeInFile=find(T>=timerange(1) & T<timerange(2));
        if isempty(TimeInFile)==1;
            EX2D=empty_structure_return(vars);
            return
        end
    end
    
    if nc_isvar(filename,'exact_time')==1;
        eT=nc_varget(filename,'exact_time');
    else
        eT=repmat(T,1,size(longitude,2));
    end

    %% now get  requested variables
 %% get general variables
    % lat, lon and time alread downloaded above
    EX2D.x=longitude(TimeInFile,:);
    EX2D.y=latitude(TimeInFile,:);
    EX2D.t=eT(TimeInFile,:);
    EX2D.tvec=T(TimeInFile);

    % use var_EX2Dist to get only variables found in the file

    if ischar(vars)==1;  % if only one variable was given
        temp_variable=nc_varget(filename,vars);
        EX2D.(vars)=temp_variable(TimeInFile,:);
        clear temp_variable
    elseif iscell(vars)==1 % if a cell of variables was given
        for K=1:length(var_exist);
            if var_exist(K)==1
                temp_variable=nc_varget(filename,vars{K});
                if min(size(temp_variable))>1
                EX2D.(vars{K})=temp_variable(TimeInFile,:);
                clear temp_variable
                end
            else
                EX2D.(vars{K})=nan(size(EX2D.t));
                %warning([vars{K} 'does not EX2Dist in' filename])
            end
        end
    end % end get requested variables


   
   

    %
    if nc_isvar(filename,'pressure')==1
        temp_variable=nc_varget(filename,'pressure');
         EX2D.zvec=temp_variable;
         EX2D.z=ones(size(EX2D.t));
         for Z=1:size(EX2D.t,1);EX2D.z(Z,:)=temp_variable;end 
       
    elseif nc_isvar(filename,'pressure')==0
        if nc_isvar(filename,'depth')==1  % change depth to pressure
           temp_variable=nc_varget(filename,'depth');
            pres = sw_pres(temp_variable,EX2D.y(1,:));
             EX2D.z=ones(size(EX2D.t));
         for Z=1:size(EX2D.t,1);EX2D.z(Z,:)=pres;end 
       
        else

            EX2D=empty_structure_return(vars);

        end
    end



end
end
%%
function [EX2D] = empty_structure_return(vars)

EX2D.x=nan;
EX2D.y=nan;
EX2D.z=nan;
EX2D.t=nan;
EX2D.tvec=nan;
EX2D.zvec=nan;

if ischar(vars)==1;
    EX2D.(vars)=nan;
elseif iscell(vars)==1;
    for M=1:length(vars);
        EX2D.(vars{M})=nan;
    end
end

end
%%






