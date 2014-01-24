function ci = obs_interpCasts(casts,zi);
%{
works on the output from obs_separateCasts, interpolating all variables to the
depth zi. zi can be a vector.
%}
%neil banas, apr 2012

for i=1:length(casts)
	c = casts(i);
	fields = fieldnames(c);
	[zu,ui] = unique(c.z);
	for fi = 1:length(fields)
		if length(zu) < 2 & length(unique(c.(fields{fi}))) == 1
			ci.(fields{fi})(i) = c.(fields{fi})(1);
		else
			try
				f = find(~isnan(c.(fields{fi})(ui)));
				ci.(fields{fi})(i) = interp1(zu(f),c.(fields{fi})(ui(f)),zi);
			catch
				ci.(fields{fi})(i) = nan;
			end
		end
	end
end