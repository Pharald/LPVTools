function b = isrow(M)
% ISROW   True if input is a row vector PMAT object
%
% ISROW(M) returns 1 (true) if the PMAT M is a row vector, i.e. 
% S=SIZE(M) has S(1) == 1.
%
% See also: isrow,  isscalar, isvector, iscolumn.

b = size(M,1)==1;

