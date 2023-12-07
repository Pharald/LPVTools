function [K,CL,GAM,INFO] = ncfsyn(P,varargin)
% NCFSYN Pointwise H-infinity normalized coprime factor synthesis for a PSS.
%
% [K,CL,GAM,INFO]=ncfsyn(G,W1,W2) performs an H-infinity normalized coprime 
% factor controller synthesis at each point in the domain of G. See
% LTI/NCFSYN for details.
%
% See also: ncfsyn, hinfsyn, mixsyn, h2syn, loopsyn.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(1, 5)
nargoutchk(0, 4)

% Array dimensions currently not allowed
niv = P.Domain.NumIV;
nad = numel(size(P))-2;
if nad>0
    error('NCFSYN does not allow PSS with array dimensions.');
end

% Call NCFSYN at each point in the domain
szP = privatesize(P);
PData = P.DataPrivate;
Domain = P.DomainPrivate;

% Loop over domain
npts = prod(szP(3:end));
K = ss( zeros([szP(2) szP(1) szP(3:end)]) );
CL = ss( zeros([szP(1)+szP(2) szP(1)+szP(2) szP(3:end)]) );
GAM = zeros([1 1 szP(3:end)]);
for i=1:npts
    [K(:,:,i),CL(:,:,i),GAM(:,:,i),INFO(i)] = ...
        ncfsyn(PData(:,:,i),varargin{:});
end

% Pack up data
CL = pss(CL, Domain);
GAM = pmat(GAM, Domain);    
K = pss(K, Domain);
INFO = pstruct( reshape(INFO,[1 1 szP(3:end)]) , Domain);

