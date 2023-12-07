function out = mldivide(A,B)
% MLDIVIDE  Left matrix divide for PMAT objects
%
% MLDIVIDE(A,B) is the result of A\B at each point in the combined
% domains of A and B.
%
% See also: mldivide, ldivide, mrdivide, rdivide, inv, mtimes.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mldivide');




% ----- OLD CODE

% if isa(a,'pmat') && isa(b,'pmat')
%     sza = size(a);
%     szb = size(b);
%     if sza(1)==szb(1)
%         out = binop(a,b,'mldivide');
%     elseif sza(1)==1 && sza(2)==1
%         aext = diag( repmat(a,[1 szb(1)]) );
%         out = mldivide(aext,b);        
%         %out = mldivide(diag(a*ones(1,szb(1))),b);
%     else
%         error('Invalid I/O Dimensions');
%     end
% else
%     out = mldivide(pmat(a),pmat(b));
% end
