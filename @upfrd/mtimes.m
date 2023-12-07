function out = mtimes(A,B)
% MTIMES  Multiply for UPFRD objects
%
% MTIMES(A,B) is the result of A*B at each point in the combined
% domains of A and B. For UPFRD objects A and B, this is equivalent to
% connecting A and B in series at each point in the combined domains.
%
% See also:  mtimes, mldivide, mrdivide, inv.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mtimes');

