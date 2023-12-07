function out = mtimes(A,B)
% MTIMES  Multiply for PSS objects
%
% MTIMES(A,B) is the result of A*B at each point in the combined
% domains of A and B. For PSS objects A and B, this is equivalent to
% connecting A and B in series at each point in the combined domains.
%
% See also:  mtimes, series, mldivide, mrdivide, inv.

% Check # of input arguments
%error(nargchk(2, 2, nargin, 'struct')) %COMPAT FIX 2016a -JT22.06.2016
narginchk(2, 2)
out = binop(A,B,'mtimes');
