function EXM = obs_merge(varargin)
%                        (EX1,EX2...type)
%                       ({EX1,EX2..},type)
%{
data can be a cell of structures or input multiple structures,
if no type is called type = 'any'
merges structures in EX1,EX2,EX3...ect
if EX1 has fields salinity, x,y,t... and EX2 has fields temperature,x,y,t
resultant EXM has all fields from input sturctures, ie.. salinity, temperature, x,y..
with nans added as necessary.

currently only works with call of type='any'
%}
%C. Bassin amd N. Banas Jan 2010
%Edited to delete .fileid from merged files C. Bassin Dec 2010 


if ischar(varargin{end})~=1    
type='any'; %default
     if iscell(varargin{1}) % works with cell of structures
        interiordata=varargin{1};
    else
        interiordata=varargin;
    end
else
type=varargin{end};
    if iscell(varargin{1})
        interiordata=varargin{1};
    else
        interiordata=varargin(1:end-1);
    end
end




switch type
    case 'any'
        F=[];
        for i=1:length(interiordata);
            F=[F;fieldnames(interiordata{i})];
        end
            AllF=unique(F);
       
            
            for k=1:length(AllF);
                 w=[];
                 for j=1:length(interiordata);
                        sizeoffield=length(interiordata{j}.x);
                     if isfield(interiordata{j},AllF{k});
                         GF=getfield(interiordata{j},AllF{k});
                         w=[w;GF];
                     else
                        w=[w;nan(sizeoffield,1)]; 
                     end
                 end
                 EXM.(AllF{k})=w;
                       
            end
            
            
     
   
    
    
    
    
    case 'all'
        warning('type=all.. currently doesnt work, resetting type = any');
        varargin{end}='any';
        EXM = obs_merge(varargin{:});
end

%get rid of fileid as it doesnt change when merged and cause incorrect
%information
if ~isfield(EXM,'filem')
    if isfield(EXM,'fileid')
    EXM = rmfield(EXM, 'fileid');
    end
end



end


