function out = lft(A,B,varargin)
% LFT  Generalized feedback interconnection of PSS objects.
% 
% M = LFT(M1,M2,NU,NY) forms the feedback interconnection of M1 above M2.
% NU specifies the number of outputs of M2 that are fed into M1.
% NY specifies the number of outputs of M1 that are fed into M2.
%
% See also: lft, feedback.
    
% lft(A,B,dim1,dim2) or lft(A,B) or lft(A,B,'name')
out = binop(A,B,'lft',varargin{:});

