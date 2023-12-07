function [V,I] = min(X,Y,DIM)
% MIN  Smallest component of PMAT objects
%
% For vectors, MIN(X) is the smallest element of X at each point in the
% domain of X. For matrices, MIN(X) is a row vector containing the minimum
% element from each column at each point in the domain of X.
%
% [Y,I]=MIN(X) returns the indices I of the minimum value at each point
% in the domain of X.
%
% MIN(X,Y) returns a PMAT with the smallest elements taken from X or
% Y at each point in the combined domains of X and Y.
%
% [Y,I] = MIN(X,[],DIM) operates along the dimension DIM. DIM can be
% 1 (row) or 2 (column).  See LPVMIN to perform a min over the indepedent
% variable dimensions.
%
% See also: min, lpvmin, max, median, mean.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 3)
nargoutchk(0, 2)

if nin==1
   DIM = 1;
   Y = [];
end

if nin==2
   % Define X and Y on a common domain
   [Xext,Yext] = domunion(X,Y);
   if nout==2
       % This should error out
       [VData,IData]=min(Xext.Data,Yext.Data);
   else
       VData=min(Xext.Data,Yext.Data);
   end
   V = pmat(VData,Xext.Domain);
elseif nin==1 || nin==3
   Data = X.Data;
   if DIM==1 || DIM==2   
      % Y should be empty
      [VData,IData]=min(Data,Y,DIM);
   else
      % Y should be empty (AD are at the end of Data)
      npiv = X.Domain.NumIV;
      [VData,IData]=min(Data,Y,DIM+npiv);
   end
   V = pmat(VData,X.Domain);
   I = pmat(IData,X.Domain);   
end

