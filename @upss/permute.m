function out = permute(m, index)
% PERMUTE  Permute row/column dimensions of a UPSS object.
%
% B = PERMUTE(A,ORDER) arranges the row/column dimensions of A to be in
% the order specified by ORDER.  B = A if ORDER is [1 2] and B = A.' if 
% ORDER is [2 1].  It is not possible to permute between the row/column
% dimensions of A and the extra dimensions associated with the independent
% variables. 
%
% See also: permute, transpose, ctranspose.

% Number of permute dimensions
np = length(index);

% Reorder as [AD IV]
niv = m.Domain.NumIV;
nad = numel(size(m))-2;
nad = max(nad,np);
Data = permute(m.Data,[(niv+1:niv+nad) (1:niv)]);

% Permute on AD
Data = permute(Data,[index nad+1:nad+niv]);

% Reorder as [IV AD]
Data = permute(Data,[(nad+1:nad+niv)  (1:nad)]);

% Pack result 
out = upss(Data,m.Domain);
