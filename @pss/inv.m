function out = inv(A)
% INV  Inverse of a PSS objects
%
% INV(A) computes the inverse of A at each point in the domain.
%
% See also: mldivide, mrdivide, mtimes.

% Check # of input arguments
narginchk(1, 1)

out = pss(inv(A.Data),A.Domain);