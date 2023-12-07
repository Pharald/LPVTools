
function [K,CL,GAM,INFO] = h2syn(P,NMEAS,NCON)
% H2SYN Pointwise H2 controller synthesis for PSS
%
% [K,CL,GAM,INFO] = h2syn(P,NMEAS,NCON) performs an H2 design at each point 
% in the domain of G. See LTI/H2SYN for details.
%
% See also: h2syn, hinfsyn, mixsyn, ncfsyn, loopsyn.

% Input error checking
nin = nargin;
nout = nargout;
narginchk(3, 3)
nargoutchk(0, 4)

% Array dimensions currently not allowed
niv = P.Domain.NumIV;
nad = numel(size(P))-2;
if nad>0
    error('H2SYN does not allow PSS with array dimensions.');
end

% Get Data
szP = privatesize(P);
PData = P.DataPrivate;
Domain = P.DomainPrivate;

% Loop over domain
npts = prod(szP(3:end));
K = ss( zeros([NCON NMEAS szP(3:end)]) );
CL = ss( zeros([szP(1)-NMEAS szP(2)-NCON szP(3:end)]) );
GAM = zeros([1 1 szP(3:end)]);
for i=1:npts
    [Ktmp,CLtmp,GAMtmp,INFOtmp]=h2syn(PData(:,:,i),NMEAS,NCON);
    if isempty(Ktmp) || isempty(GAMtmp)
       error(['Controller Synthesis failed at grid point #' num2str(i)])
    end
    K(:,:,i) = Ktmp;
    CL(:,:,i)=CLtmp;
    GAM(:,:,i)=GAMtmp;
    INFO(i) = INFOtmp;
end

% Pack up data
CL = pss(CL, Domain);
GAM = pmat(GAM, Domain);    
K = pss(K, Domain);
INFO = pstruct( reshape(INFO,[1 1 szP(3:end)]) , Domain);


