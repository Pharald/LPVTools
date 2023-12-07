function out = rcond(M)
% RCOND   Reciprocal condition estimator for PMAT objects.
%
% RCOND(M) is the reciprocal condition estimate of M at each point in 
% the domain of M. 
%
% See also: rcond, cond, norm.

szm = privatesize(M);
out = zeros([1 1 szm(3:end)]);
Data = M.DataPrivate;
for i=1:prod(szm(3:end))
    out(1,1,i) = rcond(Data(:,:,i));
end
out = pmat(out,M.DomainPrivate);
