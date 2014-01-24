function cast = obs_identifyCasts(EX)

XYT = [EX.t, EX.x, EX.y];
if sum(isnan(XYT))==3
    cast=nan;
    return
end

Unixyt = unique(XYT,'rows');
cast = nan(size(XYT,1),1);
C=1;
for k=1:size(Unixyt,1);
    cast(XYT(:,1)==Unixyt(k,1)  & XYT(:,2)==Unixyt(k,2)&  XYT(:,3)==Unixyt(k,3))=C;
    C=C+1;
end
    
end
    
