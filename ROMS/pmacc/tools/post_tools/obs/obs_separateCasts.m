function Casts = obs_separateCasts(EX)
%{
makes a new vector of structures such that EX is split into a structure of cast data
cast(1) has fields x,y,z,t,varname1,varname2...etc for one cast
cast(2) has same fields but possibly a different length
cast(3) ...etc...
%}
%C. Bassin and N. Banas Jan 2010

if length(EX.t)==0
	Casts = [];
	return;
end

if exist('EX.cast')~=1
    EX.cast = obs_identifyCasts(EX);
end

    S=unique([EX.cast]);
    
    vars=fieldnames(EX);
    for l=1:length(vars)
    Casts.(vars{l})=[];
    end
    Casts(max(S)).tempfield=nan;% make structure big
    
    
    for k=1:length(S);
        place=find(EX.cast==S(k));
        for m=1:length(vars)
        Casts(k).(vars{m})=EX.(vars{m})(place);
        end
    end
    
    Casts = rmfield(Casts, 'tempfield');
end
        