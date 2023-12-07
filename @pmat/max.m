function [V,I] = max(X,Y,DIM)
% MAX  Largest component of PMAT objects
%
% For vectors, MAX(X) is the largest element of X at each point in the
% domain of X. For matrices, MAX(X) is a row vector containing the maximum
% element from each column at each point in the domain of X.
%
% [Y,I]=MAX(X) returns the indices I of the maximum value at each point
% in the domain of X.
%
% MAX(X,Y) returns a PMAT with the largest elements taken from X or
% Y at each point in the combined domains of X and Y.
%
% [Y,I] = MAX(X,[],DIM) operates along the DIMension DIM. DIM can be
% 1 (row) or 2 (column).  See LPVMAX to perform a max over the indepedent
% variable DIMensions.
%
% See also: max, lpvmax, min, lpvmin, median, mean.

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
       [VData,IData]=max(Xext.Data,Yext.Data);
   else
       VData=max(Xext.Data,Yext.Data);
   end
   V = pmat(VData,Xext.Domain);
elseif nin==1 || nin==3
   Data = X.Data;
   if DIM==1 || DIM==2   
      % Y should be empty
      [VData,IData]=max(Data,Y,DIM);
   else
      % Y should be empty (AD are at the end of Data)
      npiv = X.Domain.NumIV;
      [VData,IData]=max(Data,Y,DIM+npiv);
   end
   V = pmat(VData,X.Domain);
   I = pmat(IData,X.Domain);   
end

