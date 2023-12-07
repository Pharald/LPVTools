function out = permute(m, index)
% PERMUTE  Permute row/column dimensions of a PSS object.
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
out = pss(Data,m.Domain);



% if isequal(sort(index),1:max(index)) && max(index)>=ndims(m)
%    [pvi,ai] = upidx(m.DomainPrivate);
%    domPV =  rgrid(m.DomainPrivate.IVName(pvi),m.DomainPrivate.IVData(pvi));
%    newindex = ai(index);
% 
%    szm = size(m); % SIZE only has to do with Array dims
%    szmnew = szm([1 2 2+index]);
%    domA = rgrid(szmnew(3:end));
%    
%    newindex = [newindex(:); pvi(:)]'; %1 through length(index)
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
%    fidx =  [indgoes];
%    permuteidx = newindex(fidx);
% 
%    out = pss(permute(m.DataPrivate,permuteidx),ugrid);
% end
% 
