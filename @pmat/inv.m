function out = inv(mat)
% INV   Matrix inverse for PMAT objects.
%
% INV(M) is the matrix inverse of M at each point in the domain of M.
%
% See also: inv, slash, pinv, cond.
mat = pmat(mat);     % PROMOTE TO PMAT
out = mat;
szm = privatesize(mat);
for i=1:prod(szm(3:end))
    out.DataPrivate(:,:,i) = inv(mat.DataPrivate(:,:,i));
end




% ----- OLD CODE

% out = mat;
% 
% szm = size(mat);
% if szm(1)==szm(2)
%     for i=1:prod(szm(3:end))
%         out.DataPrivate(:,:,i) = inv(mat.DataPrivate(:,:,i));
%     end
% else
%     error('PMAT must be square.')
% end
% 
