function out = permute(m, index)
% PERMUTE  Permute row/column dimensions of a PMAT object.
%
% B = PERMUTE(A,ORDER) arranges the row/column dimensions of A to be in
% the order specified by ORDER.  B = A if ORDER is [1 2] and B  =A.' if 
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

% Pack result as a PMAT
out = pmat(Data,m.Domain);




% ----- OLD CODE
%
% if isequal(sort(index),1:max(index)) && max(index)>=ndims(m)
%    [pvi,ai] = upidx(m.DomainPrivate);
%    domPV =  rgrid(m.DomainPrivate.IVName(pvi),m.DomainPrivate.IVData(pvi));
%    i1 = find(index==1);
%    i2 = find(index==2);
%    exd = find(index>=3);
%    newindex(exd) = ai(index(exd)-2) + 2;
% 
%    szm = size(m); % SIZE only has to do with Array dims
%    szmnew = szm(index);
%    domA = rgrid(szmnew(3:end));
%    
%    newindex([i1 i2]) = [1 2];
%    newindex = [newindex pvi'+2];
%    
%    [ugrid,in1,in2] = domunion(domPV,domA);
%    ndomPV = numel(domPV.IVName);
%    ndomA = numel(domA.IVName);
%    [~,indPV] = sort(in1);
%    [~,indA] = sort(in2);
%    PVgoes = indPV(1:ndomPV);
%    Agoes = indA(1:ndomA);
% 
%    indgoes([PVgoes]) = ndomA+(1:ndomPV);
%    indgoes([Agoes]) = 1:ndomA;
%    fidx =  [1 2 indgoes+2];
%    permuteidx = newindex(fidx);
% 
%    out = pmat(permute(m.DataPrivate,permuteidx),ugrid);
% end

