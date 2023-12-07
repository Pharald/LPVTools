function b = ctranspose(a)
% CTRANSPOSE  Complex conjugate transpose for UPSS objects.
%
% B=CTRANSPOSE(A) returns the complex conjugate transpose of A at
% each point in the domain of A.
%
% For continuous-time UPSS, if A has the transfer function H(s) at a point 
% in the domain then B is the system with transfer function H(-s).'
%
% For discrete-time UPSS, if A has the transfer function H(z) at a point 
% in the domain then B is the system with transfer function H(1/z).'
%
% See also: ctranspose, transpose, permute.

b = a;
b.DataPrivate = b.DataPrivate';
