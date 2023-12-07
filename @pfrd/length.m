function out = length(a)
% LENGTH   Length of a PFRD object.
%
% L = LENGTH(A) returns the maximum of the row and column dimension of A.
% If A is a vector PFRD, then L is the length of the vector. 
%
% See also: length.

% TODO PJS 4/29/2011: Implement numel for PFRDs?

sza = size(a);
out = max(sza(1:2));




