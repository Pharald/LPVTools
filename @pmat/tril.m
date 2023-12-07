function b = tril(a,k)
% TRIL  Extract lower triangular part of a PMAT object.
%
% TRIL(X) is the lower triangular part of a PMAT X at each point in the
% domain of X.
%
% TRIL(X,K) is the elements on and below the K-th diagonal of a PMAT X
% at each point in the domain of X.
%
% See also: tril, triu, diag.

% Input / output error checking
nin = nargin;
narginchk(1, 2)

if nin==1
    k=0;
end
[a,k] = domunion(a,k);
sza = privatesize(a);
aData = a.DataPrivate;
kData = k.DataPrivate;

b = zeros(sza);
for i=1:prod(sza(3:end))
    b(:,:,i) = tril( aData(:,:,i) , kData(:,:,i) );
end
b = pmat(b,a.DomainPrivate);

