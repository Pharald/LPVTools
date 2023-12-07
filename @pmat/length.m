function out = length(a)
% LENGTH   Length of a PMAT object.
%
% L = LENGTH(A) returns the maximum of the row and column dimension of A.
% If A is a vector PMAT, then L is the length of the vector. 
%
% See also: length.

sza = size(a);
out = max(sza(1:2));


% Original code commented out on 10/10/12

% % TODO PJS 4/4/2011: Implement numel for PMATs?
% 
% out = [];
% error('Use NUMEL')




