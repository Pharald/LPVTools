function [K,CL,GAM,INFO] = lpvmixsyn(G,varargin)
% LPVMIXSYN Mixed-sensitivity synthesis for PSS
%
% lft(P,K): [z1] = [W1 0 ][S  -SG ][W2 0 ][w1]
%           [z2]   [0  W3][KS -KSG][0  W4][w2]
%
% CL: closed loop, feedback connection of G and K with additional performance
% inputs and outputs. inputs={r,d}, outputs={e,u,y}
%
% [K,CL,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4,...) calculates the controller K by
% minimising the induces L2 norm of lft(P,K) as depicted above 
%
% [K,CL,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4,Xb,Yb) calculates the controller K by
% minimising the induces L2 norm of lft(P,K) as depicted above. Xb and Yb
% are BASIS objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% See also: lpvsyn, mixsyn, loopsyn.


nin = nargin;
nout = nargout;
narginchk(2, inf)
nargoutchk(0, 4)

if ~isa(varargin(end),'basis')
Xb = [];
Yb = Xb;
else
Xb = varargin{5};
Yb = varargin{6};
end

% assign weights
We = varargin{1};
Wr = varargin{2};
Wu = varargin{3};
Wd = varargin{4};
nmeas = size(G,1); % e
ncon = size(G,2); % u

% all output signals are in the feeback channel
if ~((size(Wr,1) == size(We,2)) && size(Wr,1) == nmeas)
error('the output of W2 and input to W1 must have the same dimensions as each other and as the output from G')
end

if ~((size(Wd,1) == size(Wu,2)) && size(Wu,2) == ncon)
error('the output of W4 and input to W3 must have the same dimensions as each other and as the input to G')
end

% checks
if ~isa(Wd,'double')
    error('W4 must be static (type double)')
end

% generalised plant
Pgen = [We*Wr, -We*G*Wd, -We*G; zeros(ncon,ncon+nmeas), Wu; Wr, -G*Wd, -G];

% perform synthesis
[K,GAM,INFO] = lpvsyn(Pgen,nmeas,ncon,varagin{5:end});

% inputs: r, d 
% outputs: e, u, y
CL = [eye(nmeas) -G; K -K*G; G*K, G];

end

function [K,CL,GAM,INFO] = LOCALlpvmixsynstruc(G,W1,W2,W3,W4,Xb,Yb,opts)
% helper function to run the structured (observer based) mixed sensitivity
% problem (see Theis, Pfifer, "Observer‐based synthesis of linear
% parameter‐varying mixed sensitivity controllers", RNC 2020)

if ~isa(W2,'double')
    error('W2 must be static (type double)')
end

% scale plant for coprime factorization
Gcop = W2\G*W4;


end