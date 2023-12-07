function B = lpvsubs(A,Names,Values,imeth)
% LPVSUBS  Evaluate PMAT at points in the domain, and demote reusults to DOUBLE.
%
% B = LPVSUBS(A,NAME,VALUES,METHOD) evaluates a PMAT A at the domain points
% specified by the NAME/VALUES pair.  NAME is an N-by-1 cell array of
% characters. NAME must contain all of the independent variable names in
% A, namely A.IVName, but may also contain others but they will be ignored. 
% VALUES is an N-by-NPTS double array of the corresponding values.  B is a 
% double array, with row and column dimensions from A.  The 3rd dimension 
% of B is NPTS.  METHOD is an optional input that specifies the 
% interpolation method and can be: 'nearest', 'linear', 'spline', or 
% 'cubic'. The default is 'linear'.
%
% See also: lpvinterp, lpvsplit, lpvsample.


% Parse inputs
nin = nargin;
narginchk(3, 4)
if nin==3
    imeth = 'linear';
end

% Check interpolation method
if isempty(strmatch(imeth,{'linear','cubic','nearest','spline','interpnlinear'})),
    error('Choose an allowable interpolation method')
end

% Get PMAT Data
szA = size(A);
IVNameA = A.Domain.IVName;
IVDataA = A.Domain.IVData;
DIVDataA = A.Domain.DIVData;
NumIVA = A.Domain.NumIV;
LIVDataA = A.Domain.LIVData;

% Check Names
NumNames = numel(Names);
if numel(unique(Names))~=NumNames
    error('Specified IVNames cannot have duplicate entries');
end

[C,IA,IB] = intersect(IVNameA,Names);  %  C = A(IA) and C = B(IB).
if numel(C)~=NumIVA
    error('Specified IVNames must include all IVs in A');
end

% Since IA reorders IVNameA, it needs to be rearranged back to the original
% order, and IB must be adjusted accordingly.  Afterwards, newIB will line
% up the rows with the original ordering in A.
[~,idx] = sort(IA);
newIB = IB(idx);

% Check Values
szV = size(Values);
if szV(1)~=NumNames
    error('Dimension of specified Names does not match Values');
end

% Perform substitution
if length(szA)==2
    nad=1;
else
    nad = prod(szA(3:end));
end

szB = [szA(1) szA(2) nad szV(2)];
B = zeros(szB);
switch imeth
    case 'linear'
        % TODO PJS 5/16/2011: This currently errors out if Values is
        % outside the domain of A.  Should we allow extrapolation and/or
        % hold the nearest value? In the least we should add a more
        % descriptive error message
        
        % TODO PA 10/10/12 - Revisit to replace for loops and/or speed up
        len = 2+NumIVA;
        idx = cellstr(repmat(':',len,1))';       
        for j = 1:nad
            for i=1:szV(2)
                Atemp = A.Data(idx{:},j);
                B(:,:,j,i) = ndLinInterp(Atemp,NumIVA,LIVDataA,IVDataA,...
                    DIVDataA,Values(newIB,i));
            end
        end
        B = reshape(B,[szA szV(2)]);
    otherwise
        % TODO: This will require a different call
        error('Only linear interpolation is supported as of yet.');
end

