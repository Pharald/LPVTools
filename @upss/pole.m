function P = pole(S)
% POLE  Pointwise computes the poles of a UPSS object
%
%  P = POLE(SYS) returns the poles P of the UPSS SYS at each point in the
%  domain of SYS.  P is returned as a column vector PMAT containing the
%  eigenvalues of the state matrix A at each point in the domain of SYS.
%
% See also: pole, damp, zero.

szS = privatesize(S);
SData = S.DataPrivate;
nx = size(SData(:,:,1).a,1);

P = zeros([nx 1 szS(3:end)]);
npts = prod(szS(3:end));
for i=1:npts
   P(:,:,i) = pole( SData(:,:,i) );
end
P = pmat(P,S.DomainPrivate);