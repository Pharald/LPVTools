function [K,CL,GAM,INFO] = loopsyn(P,varargin)
% LOOPSYN Pointwise loopsyn controller synthesis for PSS
%
% [K,CL,GAM,INFO]=loopsyn(G,Gd,..,) performs a loopsyn controller
% synthesis at each point in the domain of G. See LTI/LOOPSYN for details.
%
% See also: loopsyn, hinfsyn, mixsyn, ncfsyn, h2syn.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(2, 3)
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
CL = ss( zeros([szP(1) szP(1) szP(3:end)]) );
GAM = zeros([1 1 szP(3:end)]);
for i=1:npts
    [K(:,:,i),CL(:,:,i),GAM(:,:,i),INFO(i)] =  ...
        loopsyn(PData(:,:,i),varargin{:});
end

% Pack up data
CL = pss(CL, Domain);
GAM = pmat(GAM, Domain);    
K = pss(K, Domain);
INFO = pstruct( reshape(INFO,[1 1 szP(3:end)]) , Domain);


