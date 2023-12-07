function out = mrdivide(A,B)
% MRDIVIDE  Right matrix divide for PMAT objects
%
% MRDIVIDE(A,B) is the result of A/B at each point in the combined
% domains of A and B.
%
% See also: mrdivide, mldivide, rdivide, ldivide, inv, mtimes.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mrdivide');




% ----- OLD CODE

% if isa(a,'pmat') && isa(b,'pmat')
%     sza = size(a);
%     szb = size(b);
%     if sza(2)==szb(2)
%         out = binop(a,b,'mrdivide');
%     elseif szb(1)==1 && szb(2)==1
%         bext = diag( repmat(b,[1 sza(2)]) );
%         out = mrdivide(a,bext);
%         %out = mrdivide(a,diag(b*ones(1,sza(2))));
%     else
%         error('Invalid I/O Dimensions');
%     end
% else
%     out = mrdivide(pmat(a),pmat(b));
% end
