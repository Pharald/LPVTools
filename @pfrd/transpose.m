function b = transpose(a)
% TRANSPOSE  Transpose for PFRD objects.
%
% B=TRANSPOSE(A) returns the non-conjugate transpose of A at
% each point in the domain of A. 
% 
% If A has the transfer function H(s) at a point in the domain then B
% is the system with transfer function H(s).'
%
% See also: transpose, ctranspose, permute.

b = a;
b.DataPrivate = b.DataPrivate.';
