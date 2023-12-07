function b = iscolumn(M)
% ISCOLUMN   True if input is a column vector PMAT object
%
% ISCOLUMN(M) returns 1 (true) if the PMAT M is a column vector, i.e. 
% S=SIZE(M) has S(2) == 1.
%
% See also: iscolumn,  isscalar, isvector, isrow.

b = size(M,2)==1; %SIZE only concerned with 2nd entry
