function [out] = fill_gap(in,interp,interp_len)

%   Fills gaps (of length < interp_len) along columns with scheme specified in 'interp'
%       [out] = fill_gap(in,interp,interp_len)
%           interp -> 'linear'
%           interp_len -> max. length over which to interpolate (optional
%           if not interpolating)

	if(~strcmpi(interp,'linear'))
        fprintf('\n ERROR: ONLY LINEAR INTERPOLATION ALLOWED.');
        return
    end
    
    s = size(in);
    nan = isnan(in);
    
    for i=1:s(2)
        start = find(diff(nan(:,i)) == 1) + 1;
        stop = find(diff(nan(:,i)) == -1);
        
        % Special Cases
        if isnan(in(1,i)) == 1 % first gap starts at the beginning of var
            start(1:end+1,1)=[1;start];
        end
        if isnan(in(s(1),i)) == 1 % last gap extends till end of var
            stop(end+1,1) = s(1);
        end
        
        gaps = stop-start+1; %[start(1)-1;start(2:end)-stop(1:end-1);s(1)-stop(end)];
        g = find(gaps <= interp_len);
        %num = zeros(length(ln),1);
        out(:,i) = in(:,i);
        
%         if(isempty(g))
%             fprintf('\n No gaps filled. Depth i = %d', i);
%             continue    
%         end
        
        nn=0;
        for j=1:length(g)
            aa = start(g(j)); bb = stop(g(j));
            if aa == 1 || bb == s(1)
                continue
            end
            out(aa:bb,i) = (interp1([aa-1, bb+1],[in(aa-1,i),in(bb+1,i)],[aa:bb]))';
            nn=nn+1;
        end  
        if nn~=0, fprintf('\n Interpolated over %d of maximum %d gaps. Depth/Column %d', nn, length(gaps),i); end
    end    
    
    fprintf('\n');