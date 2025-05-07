function [K,Gamma,Info] = lpvmixsyn(G,W1,W2,W3,W4,varargin)
% LPVMIXSYN Mixed-sensitivity synthesis for PSS
%
% lft(P,K): [z1] = [W1 0 ][S  -SG ][W2 0 ][w1]
%           [z2]   [0  W3][KS -KSG][0  W4][w2]
%
%
% [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4) calculates the controller K by
% minimising the induces L2 norm of lft(P,K) as depicted above 
%
% [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4,XB,YB) calculates the controller K by
% minimising the induces L2 norm of lft(P,K) as depicted above. Xb and Yb
% are BASIS objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4,XB,YB,ALG) allows specifying either
% a full unstructured control with alg='full' or a structured controller
% with a specific filter/state-feedback structure alg ='struc'
%
% [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,W4,XB,YB,ALG,OPT) allows the user to
% pass in a LPVSYNOPTIONS object. 
%
% See also: lpvsyn, mixsyn, loopsyn.


nin = nargin;
nout = nargout;
narginchk(5, 9)
nargoutchk(0, 3)

% Parse Inputs
if nin==5
    opt = lpvsynOptions;
    alg = 'full';
    Yb = [];
    Xb = Yb;
elseif nin==6
    Yb = [];
    Xb = Yb;
    if isa(varargin{1},'lpvsynoptions')
        opt = varargin{1};
        alg = 'full';
    elseif isa(varargin{1},'char')
        alg = varargin{1};
        opt = lpvsynOptions;
    else
        error('Expected lpvsynoptions or alg as 6th input')
    end
elseif nin == 7
    if isa(varargin{1},'basis')
        Xb = varargin{1};
        Yb = varargin{2};
        opt = lpvsynOptions;
        alg = 'full';
    else
        Yb = [];
        Xb = Yb;
        alg = varargin{1};
        opt = varargin{2};
    end
elseif nin >=8
    Xb = varargin{1};
    Yb = varargin{2};
    if nin == 8
        if isa(varargin{3},'lpvsynoptions')
            opt = varargin{3};
            alg = 'full';
        elseif isa(varargin{3},'char')
            alg = varargin{3};
            opt = lpvsynOptions;
        else
            error('Expected lpvsynoptions or alg as 8th input')
        end
    end
    if nin ==9
        alg = varargin{3};
        opt = varargin{4};
    end
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



nmeas = size(G,1); % e
ncont = size(G,2); % u

% all output signals are in the feeback channel
if ~((size(W2,1) == size(W1,2)) && size(W2,1) == nmeas)
error('the output of W2 and input to W1 must have the same dimensions as each other and as the output from G')
end

if ~((size(W4,1) == size(W3,2)) && size(W3,2) == ncont)
error('the output of W4 and input to W3 must have the same dimensions as each other and as the input to G')
end



if strcmp(alg,'struc')
    
        [K,Gamma,Info] = LOCALlpvmixsynstruc(G,W1,W2,W3,W4,Xb,Yb,opt);
        
elseif strcmp(alg,'full')
    % add performance weights    
    % TODO HP 09/12/24: based on the connect command. Need to reevaluate
    % later what is the best option to build plant interconnections
    W1 = ss(W1);
    W2 = ss(W2);
    W3 = ss(W3);
    W4 = ss(W4);

    % Name Inputs/Ouputs for connect
    W1.InputName = 'e';  W1.OutputName = 'e1';
    W2.InputName = 'd1'; W2.OutputName = 'r';
    W3.InputName = 'u';  W3.OutputName = 'e2';
    W4.InputName = 'd2'; W4.OutputName ='d';
    G.InputName ='v';    G.OutputName = 'y';
    S1 = sumblk('e = r -y');
    S2 = sumblk('v = u + d');

    P = connect(G,W1,W2,W3,W4,S1,S2,{'d1','d2','u'},{'e1','e2','e'});
    % perform synthesis
    [K,Gamma,Info] = lpvsyn(P,nmeas,ncont,Xb,Yb,opt);
else
    error('alg only supports full or struc')
end


end

function [K,Gamma,Info] = LOCALlpvmixsynstruc(P,W1,W2,W3,W4,Xb,Yb,opt)
% helper function to run the structured (observer based) mixed sensitivity
% problem (see Theis, Pfifer, "Observer‐based synthesis of linear
% parameter‐varying mixed sensitivity controllers", RNC 2020)


% Parse Inputs
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
[~,~,~,Md] = ssdata(M); % -- EB 05.25 necessary for lft case until plftss objects are fixed (of M.d is a plftss instead of pmat) 
Psf = ss(A,[L B],-C,[W2/(Md) D]);

% add output performance weights 
tmp = parallel(Psf,W3,(nperf+1:nperf+ncont),1:ncont,[],[]);
Psfweighted = blkdiag(W1,eye(ncont))*tmp;

[F,Gamma,Info] = lpvsfsyn(Psfweighted,ncont,Yb,'L2',opt);

Info.L = L;

% reconstruct controller
nx = size(A,1);
nxWu = size(Wu.a,1);
nxWe = size(We.a,1);

% the state order is x_We, x_P, x_Wu
Aaug = [We.a, -We.b*C, zeros(nxWe,nxWu);
     zeros(nx,nxWe), A, zeros(nx,nxWu);
     zeros(nxWu,nxWe), zeros(nxWu,nx), Wu.a];
Laug = [ We.b ; L/W2 ; zeros(nxWu,nmeas)];
Caug = [ zeros(nmeas,nxWe), C , zeros(nmeas,nxWu)];
Baug = [ zeros(nxWe,ncont); B ; Wu.b];


K = ss(Aaug+Laug*Caug+Baug*F,Laug,F,0);


end