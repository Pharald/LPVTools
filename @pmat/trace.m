function out = trace(mat)
% TRACE  Sum of diagonal elements for PMAT objects.
%
% TRACE(M) is the sum of the diagonal elements of a square PMAT M at 
% each point in the domain of M.
%
% See also: trace.


szm = privatesize(mat);
out = zeros([1 1 szm(3:end)]);
for i=1:prod(szm(3:end))
    out(1,1,i) = trace(mat.DataPrivate(:,:,i));
end
out = pmat(out,mat.DomainPrivate);
