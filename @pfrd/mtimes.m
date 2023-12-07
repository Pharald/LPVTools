function out = mtimes(A,B)
% MTIMES  Multiply for PFRD objects
%
% MTIMES(A,B) is the result of A*B at each point in the combined
% domains of A and B. For PFRD objects A and B, this is equivalent to
% connecting A and B in series at each point in the combined domains.
%
% See also:  mtimes, mldivide, mrdivide, inv.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'mtimes');



% ----- OLD CODE

% function out = mtimes(A,B)
% 
% if isa(A,'pss') && isa(B,'pss')
%    sa = size(A);
%    sb = size(B);
%    if sa(2)==sb(1)
%       try
%          out = binop(A,B,'mtimes');
%       catch
%          error(lasterr)
%       end
%    elseif sa(1)==1 && sa(2)==1 && sb(1)>0 && sb(2)>0
%       Aext = A;
%       for i=2:sb(1)
%          Aext = blkdiag(Aext,A);
%       end
%       out = Aext*B;
%    elseif sb(1)==1 && sb(2)==1 && sa(1)>0 && sa(2)>0
%       Bext = B;
%       for i=2:sa(2)
%         Bext = blkdiag(Bext,B);
%       end
%       out = A*Bext;
%    else
%       error('Invalid inner dimensions');
%    end
% else
%    out = pss(A)*pss(B);
% end
