function V = median(X,DIM)
% MEDIAN  Median value for PMAT objects
%
% For vectors, MEDIAN(X) is the median value of the elements in X at each 
% point in the domain of X. For matrices, MEDIAN(X) is a row vector 
% containing the median value of each column at each point in the domain 
% of X.  
%
% MEDIAN(X,DIM) takes the median along the dimension DIM of X. 
%
% See also: median, lpvmedian, max, min, mean.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 1)

if nin==1
   DIM = 1;
end

Data = X.Data;
if DIM==1 || DIM==2   
   VData=median(Data,DIM);
else
   npiv = X.Domain.NumIV;
   VData=median(Data,DIM+npiv);
end
V = pmat(VData,X.Domain);

