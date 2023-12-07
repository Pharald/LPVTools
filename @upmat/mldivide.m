function out = mldivide(A,B)
% MLDIVIDE  Left division for UPMAT objects
%
% MLDIVIDE(A,B) is the result of A\B at each point in the combined
% domains of A and B.
%
% See also: mldivide, mrdivide, inv, mtimes.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mldivide');

