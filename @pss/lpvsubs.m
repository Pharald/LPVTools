function Si = lpvsubs(m,varargin)
% LPVSUBS  Evaluate PSS at points in the domain, and demote results to SS.
%
% B = LPVSUBS(A,NAME,VALUES,METHOD) evaluates a PSS A at the domain points
% specified by the NAME/VALUES pair.  NAME is an N-by-1 cell array of
% characters. NAME must contain all of the independent variable names in
% A, namely A.IVName, but may also contain others. VALUES is an N-by-NPTS
% double array of the corresponding values.  B is an SS array, with row
% and column dimensions from A.  The 3rd dimension of B is NPTS. METHOD is 
% an optional input that specifies the interpolation method and can be: 
% 'nearest', 'linear', 'spline', or 'cubic'. The default is 'linear'.
%
% See also: lpvinterp, lpvsplit, lpvsample.



% Parse inputs
nin = nargin;
narginchk(3, 4)

% Get state-space data
Data = m.DataPrivate;
[A,B,C,D]=ssdata(Data);
ns = size(A,1);
Domain = m.DomainPrivate;

% Interpolate stacked state-space matrices using PMAT/LPVINTERP
M = [pmat(A,Domain) pmat(B,Domain); pmat(C,Domain) pmat(D,Domain)];
Mi = lpvsubs(M,varargin{:});
Ai = Mi(1:ns,1:ns,:);
Bi = Mi(1:ns,ns+1:end,:);
Ci = Mi(ns+1:end,1:ns,:);
Di = Mi(ns+1:end,ns+1:end,:);

% Pack up data for output (currently delays are not handled)
% SS Syntax inherits the following generic properties from Data:
%   InputDelay, InputGroup, InputName, Notes, OutputDelay,
%   OutputGroup, OutputName, Ts, UserData
% The following properties are not documented as 'generic' but
% also seem to be set in the SS syntax below:
%   Name, TimeUnit, InputUnit, OutputUnit
Si = ss(Ai,Bi,Ci,Di,Data);

% Set non-generic properties
% TODO PJS 5/1/2011: Are there other properties to be set?
Si = set(Si,'Scaled',Data.Scaled, ...
    'StateName',Data.StateName,'StateUnit',Data.StateUnit);



% return
% 
% % Get PSS Data
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
% [Ass,Bss,Css,Dss]=ssdata(A.Data);
% 
% 
% szB = [szA(1:2) szV(2)];
% B = ss(zeros(szB));
% switch imeth
%     case 'linear'
%         % TODO PJS 5/16/2011: This currently errors out if Values is
%         % outside the domain of A.  Should we allow extrapolation and/or
%         % hold the nearest value? In the least we should add a more
%         % descriptive error message
%         for i=1:szV(2)
%             B(:,:,i) = ndLinInterp(A.Data,NumIVA,LIVDataA,IVDataA,...
%                 DIVDataA,Values(newIB,i));
%         end
%     otherwise
%         % TODO: This will require a different call
%         error('Other interpolation methods not supported yet.');
% end
% 
