function out = mrdivide(A,B)
% MRDIVIDE  Right division for PFRD objects
%
% MRDIVIDE(A,B) is the result of A/B at each point in the combined
% domains of A and B.
%
% See also: mrdivide, mldivide, mtimes.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mrdivide');
