function [K,CL,GAM,INFO] = mixsyn(P,varargin)
% MIXSYN Pointwise Mixed-sensitivity synthesis for PSS
%
% [K,CL,GAM,INFO]=mixsyn(G,W1,W2,W3,...) performs an H-infty mixed
% sensitivity synthesis at each point in the domain of G. See LTI/MIXSYN for 
% details.
%
% See also: mixsyn, hinfsyn, ncfsyn, h2syn, loopsyn.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, 4)

% Array dimensions currently not allowed
niv = P.Domain.NumIV;
nad = numel(size(P))-2;
if nad>0
    error('MIXSYN does not allow PSS with array dimensions.');
end

% Get Data
szP = privatesize(P);
PData = P.DataPrivate;
Domain = P.DomainPrivate;

% Loop over domain
npts = prod(szP(3:end));
K = ss( zeros([szP(2) szP(1) szP(3:end)]) );
GAM = zeros([1 1 szP(3:end)]);
for i=1:npts
    [K(:,:,i),CLi,GAM(:,:,i),INFO(i)] =  ...
        mixsyn(PData(:,:,i),varargin{:});
    if i==1
        CL = ss( zeros( [size(CLi) szP(3:end)] ) );
    end
    CL(:,:,i) = CLi;
end

% Pack up data
CL = pss(CL, Domain);
GAM = pmat(GAM, Domain);    
K = pss(K, Domain);
INFO = pstruct( reshape(INFO,[1 1 szP(3:end)]) , Domain);


