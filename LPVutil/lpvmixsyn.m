function [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,varargin)
% LPVMIXSYN  Parameter-varying mixed-sensitivity synthesis for PSS.
% 
% [K,GAM,INFO]=LPVMIXSYN(G,W1,W2,W3) synthesizes the parameter-varying
% mixed-sensitivity controller K, which minimizes the induced L2 norm of  
% W1*S, W2*K*S and W3*T, where S = inv(I+G*K), T = G*K*inv(I+G*K), and 
% W1, W2 and W3 are stable PSS, PMAT or DOUBLE weights of appropriate size. 
% GAM is the induced L2 norm acheived by K. INFO is a structure containing 
% data from the Linear Matrix Inequalities that are solved to obtain K. 
% A call to LPVMIXSYN without a BASIS function argument generates a 
% controller assuming no bounds on the parameter rate of variation. 
%
% [K,GAM,INFO]=LPVMIXSYN(G,W1,W2,W3,Xb,Yb) computes the rate-bounded 
% mixed-sensitivity controller K, where the rate-bounds of the independent 
% variables of G, W1, W2 and W3 are incuded in the synthesis conditions. 
% Xb and Yb are BASIS objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% [K,GAM,INFO]=LPVMIXSYN(G,...,OPT) allows the user to pass in a 
% LPVSYNOPTIONS object. 
%  
% See also: lpvsynOptions, lpvsyn, lpvsfsyn, lpvestsyn,  lpvncfsyn, lpvloopshape.

% Parse inputs
nin = nargin;
if nin==1
    W1=[];
    W2=[];
    W3=[];
    varargin = cell(0,0);
elseif nin==2
    W2=[];
    W3=[];
    varargin = cell(0,0);
elseif nin==3
    W3=[];    
    varargin = cell(0,0);
elseif nin==4
    varargin = cell(0,0);         
end

% Default weights
[ny,nu] = size(G);
Iy = eye(ny);
Iu = eye(nu);
if isempty(W1)
    W1 = ss(zeros(0,ny));
end
if isempty(W2)
    W2 = ss(zeros(0,nu));
end
if isempty(W3)
    W3 = ss(zeros(0,ny));
end

if isscalar(W1)
    W1 = W1*Iy;
end
if isscalar(W2)
    W2 = W2*Iu;
end
if isscalar(W3)
    W3 = W3*Iy;
end


% Build generalized interconnection structure
col1 = [Iy;zeros(nu+ny,ny);Iy];
col2 = [zeros(ny,nu);Iu;zeros(2*ny,nu)]+[-Iy;zeros(nu,ny);Iy;-Iy]*G;
allw = blkdiag(W1,W2,W3,Iy);
Pic = allw*[col1,col2];

% Synthesize controller Ks for shaped plant
[K,GAM,INFO] = lpvsyn(Pic,ny,nu,varargin{:}); 

