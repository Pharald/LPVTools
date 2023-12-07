function out = det(mat)
% DET  Determinant for PMAT objects.
%
% DET(M) is the determinant of a square PMAT M at each point in the 
% domain of M.
%
% See also: det, cond.


szm = privatesize(mat);
out = zeros([1 1 szm(3:end)]);
for i=1:prod(szm(3:end))
    out(1,1,i) = det(mat.DataPrivate(:,:,i));
end
out = pmat(out,mat.DomainPrivate);



% ----- OLD CODE

% szm = size(mat);
% if szm(1)==szm(2)
%     out = zeros(szm);
%     for i=1:prod(szm(3:end))
%         out(:,:,i) = det(mat.DataPrivate(:,:,i));
%     end
%     out = pmat(out,mat.DomainPrivate);
% else
%     error('Matrix must be square');
% end
