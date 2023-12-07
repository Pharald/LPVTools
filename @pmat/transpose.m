function b = transpose(a)
% TRANSPOSE  Non-conjugate transpose for PMAT objects.
%
% TRANSPOSE(A) returns the non-conjugate transpose of A at
% each point in the domain of A.
%
% See also: transpose, ctranspose, permute.

sza = privatesize(a);
b = a;
b.DataPrivate = permute(a.DataPrivate,[2 1 3:length(sza)]); 
