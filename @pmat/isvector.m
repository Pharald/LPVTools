function b = isvector(M)
% ISVECTOR   True if input is a vector PMAT object
%
% ISVECTOR(M) returns 1 (true) if the PMAT M is either a row or column
% vector, i.e. S=SIZE(M) has S(1) == 1 or S(2) == 1.
%
% See also: isvector, isscalar, isrow, iscolumn.

szM = size(M);
b = (szM(1)==1 || szM(2)==1);
