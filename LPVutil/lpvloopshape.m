function [K,Gamma,Info] = lpvloopshape(P,Wi,Wo,varargin)
% LPVLOOPSHAPE  Parameter-varying loop-shaping synthesis for PSS.
% 
% 
% [K,GAM,INFO]=LPVLOOPSHAPE(G,Wi,Wo) synthesizes the parameter-varying
% controller K, which minimizes the induced L2 norm for the shaped 
% plant Gs=Wo*G*Wi. GAM is the induced L2 norm of the closed-loop system 
% from [d1,d2] to [e1,e2].
%
%                       ^ e1                          ^ e2
%            ____       |    ____    ____     ____    |
%     +---->| Ks |-- + ---->| Wi |--| G  |-->| Wo |------->+ ----- 
%   - |      ----    ^       ----    ----     ----         ^     |
%     |              |                                     |     | 
%     |           d1 |                                  d2 |     |
%     |__________________________________________________________|
%
% INFO is a structure containing data from the Linear Matrix Inequalities 
% that are solved to obtain K. A call to LPVLOOPSHAPE without a BASIS 
% function argument generates a controller assuming no bounds on the 
% parameter rate of variation. 
%
% [K,GAM,INFO]=LPVLOOPSHAPE(G,Wi,Wo,Xb,Yb) computes the rate-bounded 
% controller K, where the rate-bounds of the independent variables of the 
% shaped pland Gs are incuded in the synthesis conditions. Xb and Yb are 
% BASIS objects, which describe the assumed parameter dependence of the 
% lyapunov matrices used in solving for K.
%
% [K,CL,GAM,INFO]=LPVLOOPSHAPE(G,...,OPT) allows the user to pass in a 
% LPVSYNOPTIONS object. 
%
% See also: lpvsynOptions, lpvsfsyn, lpvestsyn, lpvsyn, lpvncfsyn, lpvmixsyn.


% Parse Inputs
nin = nargin;
narginchk(1, 6)
nout = nargout;
if nin==1
    Wi = [];
    Wo = [];
elseif nin==2
    Wo = [];
end

% Default weights
[ny,nu] = size(P);
Iy = eye(ny);
Iu = eye(nu);
if isempty(Wo)
    Wo = ss(zeros(0,ny));
end
if isempty(Wi)
    Wi = ss(zeros(0,nu));
end

if isscalar(Wo)
    Wo = Wo*Iy;
end
if isscalar(Wi)
    Wi = Wi*Iu;
end


% Form design interconnection
%  [e1]     [ I   0   I ] [d1]
%  [e2]  =  [ Pn  0  Pn ] [d2]
%  [y ]     [-Pn -I -Pn ] [u]
Pnew = Wo*P*Wi;
szPnew = size(Pnew);
ni = szPnew(1);
no = szPnew(2);

Pic = [eye(ni) zeros(ni,no) eye(ni); zeros(no,2*ni+no); ...
   zeros(no,ni) -eye(no) zeros(no,ni)] + ...
   [zeros(ni,no); eye(no); -eye(no)]*Pnew*[eye(ni) zeros(ni,no) eye(ni)];

[K0,Gamma,Info] = lpvsyn(Pic,no,ni,varargin{:});
Info.K0 = K0;
K = Wi*K0*Wo;



