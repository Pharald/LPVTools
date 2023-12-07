function out = plus(A,B)
% PLUS  Plus for PFRD objects
%
% PLUS(A,B) is the result of A+B at each point in the combined
% domains of A and B. For PFRD objects A and B, this is equivalent to 
% connecting A and B in parallel.
%
% See also: plus, parallel.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'plus');



% ----- OLD CODE


% function out = plus(A,B)
% 
% if isa(A,'pss') && isa(B,'pss')
%    sa = size(A);
%    sb = size(B);
%    if sa(1)==1 && sa(2)==1 && (sb(1)>1 || sb(2)>1)
%       Aext = ones(sb(1),1)*A*ones(1,sb(2));
%       out = Aext+B;        
%    elseif sb(1)==1 && sb(2)==1 && (sa(1)>1 || sa(2)>1)
%       Bext = ones(sa(1),1)*B*ones(1,sa(2));
%       out = A+Bext;        
%    elseif sa(1)==sb(1) && sa(2)==sb(2)
%       try
%          out = binop(A,B,'plus');
%       catch
%          error(lasterr)
%       end
%    else
%       error('Invalid Dimensions');
%    end
% else
%    out = pss(A) + pss(B);
% end



