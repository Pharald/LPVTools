function V = mean(X,dim)
% MEAN  Average or mean value for PMAT objects
%
% For vectors, MEAN(X) is the mean value of the elements in X at each 
% point in the domain of X. For matrices, MEAN(X) is a row vector 
% containing the mean value of each column at each point in the domain 
% of X.  
%
% MEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
% See also: mean, lpvmean, max, min, median.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 1)

if nin==1
   dim = 1;
end

Data = X.Data;
if dim==1 || dim==2   
   VData=mean(Data,dim);
else
   npiv = X.Domain.NumIV;
   VData=mean(Data,dim+npiv);
end
V = pmat(VData,X.Domain);

