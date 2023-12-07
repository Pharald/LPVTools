function out = permute(m, index)
% PERMUTE  Permute row/column dimensions of a UPMAT object.
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

% % Reorder as [AD IV]
% niv = m.Domain.NumIV;
% nad = numel(size(m))-2;
% nad = max(nad,np);
% Data = permute(m.Data,[(niv+1:niv+nad) (1:niv)]);
% 
% % Permute on AD
% Data = permute(Data,[index nad+1:nad+niv]);
% 
% % Reorder as [IV AD]
% Data = permute(Data,[(nad+1:nad+niv)  (1:nad)]);
% 
% % Pack result 
% out = upmat(Data,m.Domain);

% TODO AH - COPIED THE CODE FROM PMAT, BUT NOT SURE WHAT THE DESIRED 
% FUNCTIONALITY IS SUPPOSED TO BE.
%
% Reorder as [row col AD IV]
niv = m.Domain.NumIV;
nad = numel(size(m))-2;
nad = max(nad,np-2);
Data = permute(m.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);

% Permute on row, col, and AD
Data = permute(Data,[index 3+nad:2+nad+niv]);

% Reorder as [row col IV AD]
Data = permute(Data,[1 2 (3+nad:2+nad+niv)  (3:2+nad)]);

% Pack result as a UPMAT
out = upmat(Data,m.Domain);

