function out = length(a)
% LENGTH   Length of a PSTRUCT object.
%
% L = LENGTH(A) returns the maximum of the row and column dimension of A.
% If A is a vector PSTRUCT, then L is the length of the vector. 
%
% See also: length.

sza = size(a);
out = max(sza(1:2));






