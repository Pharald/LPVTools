function b = ctranspose(a)
% CTRANSPOSE  Complex conjugate transpose for PMAT objects.
%
% CTRANSPOSE(A) returns the complex conjugate transpose of A at
% each point in the domain of A.
%
% See also: ctranspose, transpose, permute.

sza = privatesize(a);
b = a;
b.DataPrivate = conj( permute(a.DataPrivate,[2 1 3:length(sza)]) ); 

