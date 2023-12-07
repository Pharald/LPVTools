function out = mtimes(A,B)
% MTIMES  Matrix multiply for PMAT objects
%
% MTIMES(A,B) is the result of A*B at each point in the combined
% domains of A and B.
%
% See also: mtimes, times, mldivide, mrdivide, inv.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mtimes');




% ----- OLD CODE

% if nargin==2
%     if isa(A,'pmat') && isa(B,'pmat')
%         sza = size(A);
%         szb = size(B);
%         if sza(2)==szb(1)
%             out = binop(A,B,'mtimes');
%         elseif sza(1)==1 && sza(2)==1
%             Aext = diag( repmat(A,[1 szb(1)]) );
%             out = Aext*B;
%             
% %             tmpa = A;
% %             for j=2:szb(1)
% %                 tmpa = blkdiag(tmpa,A);
% %             end
% %             out = tmpa*B;
%         elseif szb(1)==1 && szb(2)==1
%             Bext = diag( repmat(B,[1 sza(2)]) );
%             out = A*Bext;
%             
% %             tmpb = B;
% %             for j=2:sza(2)
% %                 tmpb = blkdiag(tmpb,B);
% %             end
% %             out = A*tmpb;
%         else
%             error('Inner matrix dimensions must agree.')
%         end
%         %out = genredu(out);
%     else
%         out = mtimes(pmat(A),pmat(B));
%     end
% end
