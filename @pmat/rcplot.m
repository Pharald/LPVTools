function H = rcplot(m)
% RCPLOT  Plot the value of PMAT matrix elements over the domain.
% 
% RCPLOT(M) plots the value of each matrix element in a PMAT over the range
% of parameter values in the PMAT's domain. If M is a R-by-C matrix that 
% depends on N indenpendent variables, then RCPLOT produces a figure with 
% R*C sublots arranged in R rows and C columns. The (i,k) subplot shows how
% the value of the matrix element in row i and column k of M changes
% as a function of the N parameters in the domain of M. RCPLOT is 
% restricted to PMATs that depend on one or two parameters. If N=1, each
% subplot plots the value of the element on the vertical axis and the value 
% of the parameter on the horizontal axis. If N =2, each subplot is a 
% surface plot with the value of the element on the vertical axis, and the
% values of the two parameters on the remaining axis.
%
% See also: plot.

% TODO PJS 4/4/2011: Revisit

szm = privatesize(m);
niv = m.DomainPrivate.NumIV;
h = zeros(szm(1),szm(2));
for i=1:szm(1)
	for j=1:szm(2)
		S(1).type = '()';
		S(1).subs = {i,j};
		subplot(szm(1),szm(2),j+(i-1)*szm(2));
		if niv==2
		    h(i,j) = surf(subsref(m,S));
		elseif niv==1
		    h(i,j) = lpvplot(subsref(m,S));
		else
			error('Invalid number of IV arguments');
		end
	end
end
if nargout==1
   H = h;
end

