function b = ctranspose(a)
% CTRANSPOSE  Complex conjugate transpose for UPMAT objects.
%
% B=CTRANSPOSE(A) returns the complex conjugate transpose of A at
% each point in the domain of A.
%
% See also: ctranspose, transpose, permute.


b = a;
b.DataPrivate = b.DataPrivate';
