function out = mtimes(A,B)
% MTIMES  Multiply for UPMAT objects
%
% MTIMES(A,B) is the result of A*B at each point in the combined
% domains of A and B. 
%
% See also:  mtimes, mldivide, mrdivide, inv.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mtimes');

