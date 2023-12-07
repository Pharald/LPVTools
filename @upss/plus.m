function out = plus(A,B)
% PLUS  Plus for UPSS objects
%
% PLUS(A,B) is the result of A+B at each point in the combined
% domains of A and B. For UPSS objects A and B, this is equivalent to 
% connecting A and B in parallel.
%
% See also: plus, parallel.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'plus');

