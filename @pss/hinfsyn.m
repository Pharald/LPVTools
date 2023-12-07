function [K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,varargin)
% HINFSYN Pointwise H-infinity controller synthesis for PSS
% 
% [K,CL,GAM,INFO] = hinfsyn(P,NMEAS,NCON,...) performs an H-infty
% design at each point in the domain of G. See LTI/HINFSYN for details.
%
% See also: hinfsyn, mixsyn, ncfsyn, h2syn, loopsyn.

% Input error checking
%nin = nargin;
%nout = nargout;
narginchk(3, inf)
nargoutchk(0, 4)

% Array dimensions currently not allowed
nad = numel(size(P))-2;
if nad>0
    error('HINFSYN does not allow PSS with array dimensions.');
end

% Get Data
szP = privatesize(P);
PData = P.DataPrivate;
Domain = P.DomainPrivate;

% Handle PMAT varargin
PMATidx = find(cellfun('isclass',varargin,'pmat'));
for i = 1:numel(PMATidx)
   tmp = varargin{PMATidx(i)};
   varargin{PMATidx(i)} = tmp.Data;
end

% Loop over domain
npts = prod(szP(3:end));
K = ss( zeros([NCON NMEAS szP(3:end)]) );
CL = ss( zeros([szP(1)-NMEAS szP(2)-NCON szP(3:end)]) );
GAM = zeros([1 1 szP(3:end)]);
for i=1:npts
    vi = varargin;
    for j = 1:numel(PMATidx)
       vi{PMATidx(j)} = vi{PMATidx(j)}(:,:,i);
    end
    [Ktmp,CLtmp,GAMtmp,INFOtmp]=hinfsyn(PData(:,:,i),NMEAS,NCON,vi{:});
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
% INFO = pstruct( reshape(INFO,[1 1 szP(3:end)]) , Domain);
INFO = [];

