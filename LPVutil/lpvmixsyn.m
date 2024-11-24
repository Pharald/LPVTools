function [K,CL,Gamma,Info] = lpvmixsyn(G,varargin)
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


% assign weights
We = varargin{1};
Wr = varargin{2};
Wu = varargin{3};
Wd = varargin{4};
nmeas = size(G,1); % e
ncont = size(G,2); % u

% all output signals are in the feeback channel
if ~((size(Wr,1) == size(We,2)) && size(Wr,1) == nmeas)
error('the output of W2 and input to W1 must have the same dimensions as each other and as the output from G')
end

if ~((size(Wd,1) == size(Wu,2)) && size(Wu,2) == ncont)
error('the output of W4 and input to W3 must have the same dimensions as each other and as the input to G')
end


flgstruc = 1;

if flgstruc
    if isa(G,'pss')
        [K,Gamma,Info] = LOCALlpvmixsynstruc(G,varargin{1:end});
    else
        error('currently only grid based structured synthesis supported')
    end

else
    % generalised plant
    % TODO HP 21/11/24: double check this. It seems this code would not
    % generate a minimal realiziation of the generalized plant.
    Pgen = [We*Wr, -We*G*Wd, -We*G; zeros(ncont,ncont+nmeas), Wu; Wr, -G*Wd, -G];

    % perform synthesis
    [K,Gamma,Info] = lpvsyn(Pgen,nmeas,ncont,varargin{5:end});
end



% inputs: r, d 
% outputs: e, u, y
CL = 1; %[eye(nmeas) -G; K -K*G; G*K, G];

end

function [K,Gamma,Info] = LOCALlpvmixsynstruc(P,W1,W2,W3,W4,Xb,Yb,opt)
% helper function to run the structured (observer based) mixed sensitivity
% problem (see Theis, Pfifer, "Observer‐based synthesis of linear
% parameter‐varying mixed sensitivity controllers", RNC 2020)


% Parse Inputs
nin = nargin;
narginchk(5, 8)
nout = nargout;
if nin<=6
    if nin==5
        opt = lpvsynOptions;
    else
        opt = Xb;
    end
    Yb = [];
    Xb = Yb;
elseif nin==7
    opt = lpvsynOptions;
elseif nin~=8
    error('Incorrect number of input arguments.')
end

if isempty(Yb)
    Yb = basis(1,0);
elseif ~isa(Yb,'basis')
    error('Yb must be a BASIS object')
end
if isempty(Xb)
    Xb = basis(1,0);
elseif ~isa(Xb,'basis')
    error('Xb must be a BASIS object')
end

% Input weights must be static (either double or pmat)
if ~isa(W4,'double')
    error('W4 must be static (type double)')
end

if ~isa(W2,'double')
    error('W2 must be static (type double)')
end

We = ss(W1);
Wu = ss(W3);

% scale plant for coprime factorization
Pcop = W2\P*W4;

[~,M,~,L] = lpvlccf(Pcop,Xb);

[A,B,C,D] = ssdata(P);
[nmeas,ncont] = size(P); % u
nperf = size(M,2);

% build plant for state feedback synthesis
Psf = ss(A,[L B],-C,[W2/(M.d) D]);

% add output performance weights 
% TODO HP 24/11/24: the ordering of the states of the weighted plant do not
% match with the controller reconstruction. Name them and use connect or
% build it up explicitly
tmp = parallel(Psf,W3,(nperf+1:nperf+ncont),1:ncont,[],[]);
Psfweighted = blkdiag(W1,eye(ncont))*tmp;

[F,Gamma,Info] = lpvsfsyn(Psfweighted,ncont,Yb,'L2',opt);

% reconstruct controller
nx = size(A,1);
nxWu = size(Wu.a,1);
nxWe = size(We.a,1);

Aaug = [A, zeros(nx, nxWe+nxWu); 
        -We.b*C, We.a, zeros(nxWe,nxWu);
        zeros(nxWu,nx), zeros(nxWu,nxWe), Wu.a ];
Laug = [L/W2 ; We.b ; zeros(nxWu,nmeas)];
Caug = [C, zeros(nmeas,nxWe), zeros(nmeas,nxWu)];
Baug = [B ; zeros(nxWe,ncont); Wu.b];


K = ss(Aaug+Laug*Caug+Baug*F,Laug,F,0);


end