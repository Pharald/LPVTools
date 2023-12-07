function arg = end(mat,slot,nslots)
% END   Overloaded END for UPMAT objects.
%
% END(MAT,DIM,NUMINDICES) returns the index corresponding to the last
% entry along the dimension DIM in the UPMAT.  NUMINDICES is the number of
% indices used in the indexing expression.
%
% See also: end.

szm = size(mat); % Use SIZE since only working on Array dims
if slot < nslots
    arg = szm(slot);
else
    if slot ==1
        error('NUMINDICES must be  larger than 1 for uncertain objects.');
    elseif slot == 2
        arg = szm(slot);
    else
        % Handle mat(3:end) or mat(:,:,:,:,4:end) where remaining
        % dimensions are aggregated into last index.
        arg = prod(szm(slot:end));
    end
end






% % OLD CODE
% szm = size(mat);
% if nslots <=2
%     if nslots==2
%         arg = szm(slot);
%     else
%         arg = szm(1)*szm(2);
%     end
% else
%     error('NUMINDICES must be less than or equal to 2.');
% end
