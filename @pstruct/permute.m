function out = permute(m, index)
% PERMUTE  Permute row/column dimensions of a PSTRUCT object.
%
% B = PERMUTE(A,ORDER) arranges the row/column dimensions of A to be in
% the order specified by ORDER.  B=A if ORDER is [1 2] and B=A.' if 
% ORDER is [2 1].  It is not possible to permute between the row/column
% dimensions of A and the extra dimensions associated with the independent
% variables. 
%
% See also: permute, transpose, ctranspose.

% Number of permute dimensions
np = length(index);

% Reorder as [row col AD IV]
niv = m.Domain.NumIV;
nad = numel(size(m))-2;
nad = max(nad,np-2);
Data = permute(m.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);

% Permute on row, col, and AD
Data = permute(Data,[index 3+nad:2+nad+niv]);

% Reorder as [row col IV AD]
Data = permute(Data,[1 2 (3+nad:2+nad+niv)  (3:2+nad)]);

% Pack result as a PSTRUCT
out = pstruct(Data,m.Domain);


