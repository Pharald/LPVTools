function out = minus(A,B)
% MINUS  Minus for UPFRD objects
%
% MINUS(A,B) is the result of A-B at each point in the combined
% domains of A and B.
%
% See also: minus, plus.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'minus');

