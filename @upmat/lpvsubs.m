function B = lpvsubs(A,varargin)
% LPVSUBS  Evaluate UPMAT at point in the domain, and demote results to UMAT.
%
% B = LPVSUBS(A,NAME,VALUES,METHOD) evaluates a UPMAT A at the domain points
% specified by the NAME/VALUES pair.  NAME is an N-by-1 cell array of
% characters. NAME must contain all of the independent variable names in
% A, namely A.IVName, but may also contain others. VALUES is an N-by-NPTS
% double array of the corresponding values.  B is a double array, with row
% and column dimensions from A.  The 3rd dimension of B is NPTS. METHOD
% is an optional input that specifies the interpolation method and can be: 
% 'nearest', 'linear', 'spline', or 'cubic'. The default is 'linear'.
%
% See also: lpvinterp, lpvsplit, lpvsample, usubs.


% Parse inputs
% nin = nargin;
% error(nargchk(3, 4, nin, 'struct'))
% if nin==3
%     imeth = 'linear';
% end

%  Get state-space data
try
   [M,Delta] = lftdata(A.Data);
catch
   error('Cannot interpolate matricies whose uncertainty description changes');
end
pM = pmat(M,A.Domain); 
B = lft(Delta,lpvsubs(pM,varargin{:}));

% 
% % Check interpolation method
% if isempty(strmatch(imeth,{'linear','cubic','nearest','spline','interpnlinear'})),
%     error('Choose an allowable interpolation method')
% end
% 
% % Get UPMAT Data
% szA = [size(A) 1];
% IVNameA = A.DomainPrivate.IVName;
% IVDataA = A.DomainPrivate.IVData;
% DIVDataA = A.DomainPrivate.DIVData;
% NumIVA = A.DomainPrivate.NumIV;
% LIVDataA = A.DomainPrivate.LIVData;
% szB = szA;
% 
% % Check Names
% NumNames = numel(Names);
% if numel(unique(Names))~=NumNames
%     error('Specified IVNames cannot have duplicate entries');
% end
% 
% [C,IA,IB] = intersect(IVNameA,Names);  %  C = A(IA) and C = B(IB).
% if numel(C)~=NumIVA
%     error('Specified IVNames must include all IVs in A');
% end
% 
% % Since IA reorders IVNameA, it needs to be rearranged back to the original
% % order, and IB must be adjusted accordingly.  Afterwards, newIB will line
% % up the rows with the original ordering in A.
% [~,idx] = sort(IA);
% newIB = IB(idx);
% 
% % Check Values
% szV = size(Values);
% if szV(1)~=NumNames
%     error('Dimension of specified Names does not match Values');
% end
% 
% % Perform substitution
% szB = [szA(1:2) szV(2)];
% B = zeros(szB);
% switch imeth
%     case 'linear'
%         % TODO PJS 5/16/2011: This currently errors out if Values is
%         % outside the domain of A.  Should we allow extrapolation and/or
%         % hold the nearest value? In the least we should add a more
%         % descriptive error message
%         for i=1:szV(2)
%             B(:,:,i) = ndLinInterp(A.DataPrivate,NumIVA,LIVDataA,IVDataA,...
%                 DIVDataA,Values(newIB,i));
%         end
%     otherwise
%         % TODO: This will require a different call
%         error('Other interpolation methods not supported yet.');
% end
% 
