function out = plus(A,B)
% PLUS  Plus for UPMAT objects
%
% PLUS(A,B) is the result of A+B at each point in the combined
% domains of A and B. 
%
% See also: plus.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'plus');


