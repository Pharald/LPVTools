function out = length(a)
% LENGTH   Length of a UPMAT object.
%
% L = LENGTH(A) returns the maximum of the row and column dimension of A.
% If A is a vector UPMAT, then L is the length of the vector. 
%
% See also: length.

% TODO PJS 4/29/2011: Implement numel for UPMATs?

sza = size(a);
out = max(sza(1:2));




