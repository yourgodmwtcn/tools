function EX = obs_omit(EX,varargin);

% function EX = obs_omit(EX,criterion
%                       ...'t',criterion
%                       ...'z','criterion
% eliminates points from all fields where some criterion is met.
% criterion can be any boolean field, e.g.
%
%    EX = obs_omit(EX, EX.salinity < 0 | isnan(EX.salinity));
%
% an older syntax is also supported:
%    for criterion, just use filedname such as salinity, not full structure name EX.salinity
%        EX = obs_omit(EX, 'salinity < 0 ');
%        EX = obs_omit(EX, 'salinity == 2 ');
%    special case, getting rid of nans
%        EX = obs_omit(EX, 'salinity == nan');
%    criterion must be string
%
% for 2 dimensional data, will get rid of entire column or row associated
% with bad data,
%     EX = obs_omit(EX,'z', 'salinity == 2 '); gets rid of depth associated row
%     EX = obs_omit(EX,'t', 'salinity == 2 '); gets rid of time associated column
% C. Bassin and N. Banas Oct 2010

if length(varargin)==1
    criterion=varargin{1};
    D=1;
else
    rc=varargin{1};
    criterion=varargin{2};
    D=2;
end
    

vars = fieldnames(EX);
if ischar(criterion) % the older syntax
	if ~isempty(regexp(criterion,'nan'))
		a=regexp(criterion,'\W');
		[row,col] = eval(['find(isnan(EX.' criterion(1:a(1)-1) ')==1);']);
	else
		[row,col] = eval(['find( EX.' criterion ');']);
	end
else
	[row,col] = find(criterion);
end
 switch D
     case 1  % 1 dimensional data
        for l=1:length(vars)
            EX.(vars{l})(row)=[];
        end
     case 2  % 2 dimensional data  % need to fix for tvec and zvec
        for l=1:length(vars)
            if strcmp(rc,'t')
               if min(size(EX.(vars{l}))>1); 
                EX.(vars{l})(:,col)=[];
               end
            
            elseif strcmp(rc,'z')
                if min(size(EX.(vars{l}))>1); 
                EX.(vars{l})(row,:)=[];
                end
            
            end
        end
        
        if strcmp(rc,'t')
            EX.tvec(col)=[];
        elseif strcmp(rc,'z')   
            EX.zvec(row)=[];  
        end
            
 end

        
end