function b = triu(a,k)
% TRIU  Extract upper triangular part of a PMAT object.
%
% TRIU(X) is the upper triangular part of a PMAT X at each point in the
% domain of X.
%
% TRIU(X,K) is the elements on and above the K-th diagonal of a PMAT X
% at each point in the domain of X.
%
% See also: triu, tril, diag.

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
    b(:,:,i) = triu( aData(:,:,i) , kData(:,:,i) );
end
b = pmat(b,a.DomainPrivate);

