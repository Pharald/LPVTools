function b = transpose(a)
% TRANSPOSE  Transpose for UPMAT objects.
%
% B=TRANSPOSE(A) returns the non-conjugate transpose of A at
% each point in the domain of A. 
%
% See also: transpose, ctranspose, permute.

b = a;
b.DataPrivate = b.DataPrivate.';
