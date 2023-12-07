function out = plus(A,B)
% PLUS  Plus for PMAT objects
%
% PLUS(A,B) is the result of A+B at each point in the combined
% domains of A and B.
%
% See also: plus.

% Check # of input arguments
narginchk(2, 2)
out = binop(A,B,'plus');



% ----- OLD CODE

% if nargin==2
%     if isa(A,'pmat') && isa(B,'pmat')
%         sza = size(A);
%         szb = size(B);
%         
%         if isempty(A) || isempty(B)
%             % Handle empty cases
%             if isempty(A) && all(szb==[1 1])
%                 % empty+scalar = empty(sza)
%                 out=pmat(zeros(sza));
%             elseif isempty(B) && all(sza==[1 1])
%                 % scalar+empty = empty(szb)
%                 out=pmat(zeros(szb));
%             elseif all(sza==szb)
%                 % empty+empty = empty(sza)  [if dimensions agree]
%                 out=pmat(zeros(sza));
%             else
%                 error('Matrix dimensions must agree.');
%             end
%         elseif sza(1)==szb(1) && sza(2)==szb(2)
%             out = binop(A,B,'plus');
%         elseif sza(1)==1 && sza(2)==1 
%             %Aext = ones(szb(1),1)*A*ones(1,szb(2));
%             Aext = A;
%             niv = Aext.DomainPrivate.NumIV;
%             Aext.Data = repmat(Aext.Data,[szb(1:2) ones(1,niv)]);
%             out = Aext+B;
%         elseif szb(1)==1 && szb(2)==1 
%             %Bext = ones(sza(1),1)*B*ones(1,sza(2));
%             Bext = B;
%             niv = Bext.DomainPrivate.NumIV;
%             Bext.Data = repmat(Bext.Data,[sza(1:2) ones(1,niv)]);
%             out = A+Bext;
%         else
%             error('Matrix dimensions must agree.')
%         end
%     else
%         out = plus(pmat(A),pmat(B));
%     end
% end
