function arg = end(mat,slot,nslots)
% END   Overloaded END for PMAT objects.
%
% END(MAT,DIM,NUMINDICES) returns the index corresponding to the last
% entry along the dimension DIM in the PMAT.  NUMINDICES is the number of
% indices used in the indexing expression.
%
% See also: end.

% AH 8/29/14 - Updated to handle mat(3:end) just like a DOUBLE would.

szm = size(mat); % Use SIZE since only working on Array dims
if slot < nslots
    arg = szm(slot);
else
    % Handle mat(3:end) or mat(:,:,:,:,4:end) where remaining
    % dimensions are aggregated into last index.
    arg = prod(szm(slot:end));
end

% OLD CODE:
% szm = size(mat); % Use SIZE since only working on Array dims
% if slot==nslots && nslots>2
%    arg = prod(szm(slot:end));
% else
%    arg = szm(slot);
% end
