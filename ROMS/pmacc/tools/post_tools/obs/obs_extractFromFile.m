function EX = obs_extractFromFile(filename, vars, timerange,varargin)
%                                                ...,'section',x,y,range)
%                                               ...,'section',track,range)
%                                                 ..., 'point',x,y);
%                                                  ...,'polygon',x,y)
%                                                    ...,'all')    
%edited C. Bassin Dec 2010                                  
%%

% if (nargin~=7 ||nargin~=6|| nargin~=4)
%     error ('Not correct number of arguments')
% end

% check filename exists:
 if exist(filename, 'file')~=2
        warning(['Filename ' filename ' not found, returning empty structure']);
       EX=empty_structure_return(vars);
       return
       
 elseif exist(filename, 'file')==2
        
       %check variables exist in file
       %first check that time, lat and lon and pressure or depth exist which are necessary for
       %program to run prperly
       if nc_isvar(filename,'time')==0;
           EX=empty_structure_return(vars);
                    return
       end
       if nc_isvar(filename,'longitude')==0;
           EX=empty_structure_return(vars);
                    return
       end
       if nc_isvar(filename,'latitude')==0;
           EX=empty_structure_return(vars);
                    return
       end
       if nc_isvar(filename,'pressure')==0 
               if nc_isvar(filename,'depth')==0
                EX=empty_structure_return(vars);
                return
               end
       end
       % how many variables, check if it is a variable name or cell array, or all
        
       % if 'all'  replace vars with cell of all variables
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
            
            %no variables exist in files
                if sum(var_exist)==0
                   % warning(['no requested variables exist in the files' filename ', returning nan structure']);
                    EX=empty_structure_return(vars);
                    return
                end
                
                
             %get time array and check that time exist
                
                 T=nc_varget(filename,'time');
                              
                    if ischar(timerange)==1
                        if strcmpi(timerange,'all')==1
                            TimeInFile=1:length(T);
                        end
                    else

                        TimeInFile=find(T>=timerange(1) & T<timerange(2));
                         if isempty(TimeInFile)==1;
                             EX=empty_structure_return(vars);
                             return
                         end
                    end
                
                %{ 
                at least one variable exists, so get lat, lon, associate
                within time frame, then get variables that exist
                %}
                
                     longitude=nc_varget(filename,'longitude');
                     latitude=nc_varget(filename,'latitude');
                
                     %which type of geography?
                    switch varargin{1}
                        case 'section'  % section, trackline type
                            if isstruct(varargin{2})
                              [dist, distFrom]=trackDist(longitude(TimeInFile),latitude(TimeInFile), varargin{2});   
                               includedata=find(distFrom<=varargin{3});
                                if ~isempty(includedata)
                                 EX.dist=dist(includedata);
                                 EX.distFrom=distFrom(includedata);
                                end
                            else
                             [dist, distFrom]=trackDist(longitude(TimeInFile),latitude(TimeInFile), varargin{2},varargin{3});
                             includedata=find(distFrom<=varargin{4});
                                if ~isempty(includedata)
                                 EX.dist=dist(includedata);
                                 EX.distFrom=distFrom(includedata);
                                end
                            end
                         
                        case 'point'  % for one location  only:
                            
                               ring=make_range_ring(varargin{2}, varargin{3}, .1);%lon,;lat,,range
                              inpoly=inpolygon(longitude(TimeInFile), latitude(TimeInFile), ring(:,1), ring(:,2));
                               includedata=find(inpoly==1);
                                if ~isempty(includedata)
                                 [dist, distFrom]=trackDist(longitude(TimeInFile(includedata)),latitude(TimeInFile(includedata)), varargin{2},varargin{3});
                             EX.dist=dist;
                             EX.distFrom=distFrom;
                                end
                                                    
                     
                         case 'polygon' % polygon
                             inpoly=inpolygon(longitude(TimeInFile), latitude(TimeInFile), varargin{2}, varargin{3});
                             includedata=find(inpoly==1);
                                
                         case  'all' % all
                              includedata=1:length(TimeInFile); 
                            
                    end % end switch case type geography
                    
                                if isempty(includedata)  % if no data in area
                                    EX=empty_structure_return(vars);
                                    return
                                end
                    %% now get  requested variables
                    
                        % use var_exist to get only variables found in the file
                       
                        
                        if ischar(vars)==1;  % if only one variable was given
                                               
                                                                                     
                            temp_variable=nc_varget(filename,vars);
                            EX.(vars)=temp_variable(TimeInFile(includedata));
                            clear temp_variable
                            
                        elseif iscell(vars)==1 % if a cell of variables was given
                            for K=1:length(var_exist);
                                if var_exist(K)==1
                                    temp_variable=nc_varget(filename,vars{K});
                                    EX.(vars{K})=temp_variable(TimeInFile(includedata));
                                    clear temp_variable
                                else
                                    EX.(vars{K})=nan(length(includedata),1);
                                    %warning([vars{K} 'does not exist in' filename])
                                end
                            end
                        end % end get requested variables
                                    
                        
                     %% get general variables
                     % lat, lon and time alread downloaded above
                     EX.x=longitude(TimeInFile(includedata));
                     EX.y=latitude(TimeInFile(includedata));
                     EX.t=T(TimeInFile(includedata));
                                   
                     %
                     if nc_isvar(filename,'pressure')==1
                         temp_variable=nc_varget(filename,'pressure');
                         EX.z=temp_variable(TimeInFile(includedata));
                     elseif nc_isvar(filename,'pressure')==0
                         if nc_isvar(filename,'depth')==1  % change depth to pressure
                             temp_variable=nc_varget(filename,'depth');
                             pres = sw_pres(temp_variable(TimeInFile(includedata)),EX.y); 
                             EX.z=pres;
                         else
                            
                        EX=empty_structure_return(vars);
            
                         end
                     end
                             
                         
  
 end
end    
%%
function [EX] = empty_structure_return(vars)
    
                EX.x=nan;
                EX.y=nan;
                EX.z=nan;
                EX.t=nan;
                 
%              if ischar(vars)==1;
%                 EX.(vars)=nan;                                            
%             elseif iscell(vars)==1;
%                 for M=1:length(vars);
%                     EX.(vars{M})=nan;
%                 end
%              end
   
end       
%%
  
   
                
    
        
    
            
