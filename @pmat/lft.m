function out = lft(A,B,varargin)
% LFT  Generalized feedback interconnection of PMAT objects.
% 
% M = LFT(M1,M2,NU,NY) forms the feedback interconnection of M1 above M2
% at each point in the combined domains of M1 and M2.
% NU specifies the number of outputs of M2 that are fed into M1.
% NY specifies the number of outputs of M1 that are fed into M2.
%
% See also: lft.
    
% lft(A,B,dim1,dim2) or lft(A,B) or lft(A,B,'name')
out = binop(A,B,'lft',varargin{:});



% % ----- OLD CODE
% 
% if isa(A,'pmat') && isa(B,'pmat')
%     sza = size(A);
%     szb = size(B);
%     if nargin==2
%         if sza(1) > szb(2)
%             if szb(1) >= sza(2)
%                 error(['incorrect dimensions']);
%                 return
%             else
%                 dim1 = szb(2);
%                 dim2 = szb(1);
%             end
%         elseif sza(1) == szb(2)
%             if sza(2) == szb(1)
%                 dim1 = sza(1);
%                 dim2 = sza(2);
%             else
%                 error(['incorrect dimensions']);
%                 return
%             end
%         else
%             if szb(1) > sza(2)
%                 dim1 = sza(1);
%                 dim2 = sza(2);
%             else
%                 error(['incorrect dimensions']);
%                 return
%             end
%         end
%     elseif nargin~=4
%         error('Incorrect number of arguments');
%     end
%     if dim1 > min([sza(1) szb(2)])
%         error('incompatible dimensions')
%         return
%     elseif dim2 > min([sza(2) szb(1)])
%         error('incompatible dimensions')
%         return
%     end
%     if sza(1)==0 & sza(2)==0
%         out = B;
%     elseif szb(1)==0 & szb(2)==0
%         out = A;
%     else
%         out = binop(A,B,'starp',dim2,dim1);
%     end
% else
%     if nargin==2
%         out = lft(pmat(A),pmat(B));
%     elseif nargin==4
%         out = lft(pmat(A),pmat(B),dim1,dim2);
%     else
%         error('Incorrect number of arguments');
%     end
% end
