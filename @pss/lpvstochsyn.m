function [K,Gamma,Info] = lpvstochsyn(P,nmeas,ncont,Xb,Yb,opt)
% LPVSTOCHSYN  LPV controller synthesis for stochastic LPV systems (PSS).
%
% [K,GAMMA,INFO] = LPVSTOCHSYN(P,NMEAS,NCON) computes a parameter-varying 
% controller for the PSS P. NCON specifies the number of available control 
% inputs in P. NMEAS specifies the number of available measurements being 
% output from P. K is the parameter-dependent controller for the stochastic 
% plant P, which minimizes the stochastic LPV bound GAMMA. The stochastic 
% LPV bound is defined as the expected value of the average instantaneous 
% power of the output of P, assuming its inputs are zero mean, white-noise 
% processes with unit intensity. INFO is a struct with additional data. 
% 
% [K,GAMMA,INFO] = LPVSTOCHSYN(P,NMEAS,NCON,XB,YB) performs a rate-bounded 
% synthesis. Xb and Yb are BASIS objects which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for K.
%
% [K,GAMMA,INFO] = LPVSTOCHSYN(P,NMEAS,NCON,XB,YB,OPT) allows the user to 
% pass in a LPVSYNOPTIONS object. 
%
% The default algorithm for LPVSTOCHSYN will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% stochastic LPV bound of lft(P,K). The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% found to lie within 15% of the optimal stochastic LPV bound of lft(P,K) 
% from the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be reset to a different value in LPVSYNOPTIONS.
%
% See also: lpvsynOptions, lpvsyn, lpvsfsyn, lpvestsyn, lpvncfsyn, lpvmixsyn, lpvloopshape.




% XXX Add input parsing
alg = 'LQG';
Xb = [];
Yb = [];
opt = lpvsynOptions;

% Define the I/O sizes
szP = size(P);
ne = szP(1)-nmeas;
nd = szP(2)-ncont;

% State feedback and output estimation
L.type = '()';
L.subs = {1:ne,':'};
Psf = subsref(P,L);
[Fsf,GamSF,InfoSF] = lpvsfsyn(Psf,ncont,alg);

L.subs = {':',1:nd};
Pest = subsref(P,L);
[Lest,GamEst0,InfoEst] = lpvFCsyn(Pest,nmeas,alg);

% Controller Reconstruction:
[A,B,C,D] = ssdata(P);
D22 = D(ne+1:end,nd+1:end);
Ak = A+Lest*C(ne+1:end,:)+B(:,nd+1:end)*Fsf+Lest*D22*Fsf;
K = pss(Ak,-Lest,Fsf,zeros(ncont,nmeas));

% Estimation cost needs to be computed using transformed feedback gain
[qq,rr]=qr(D(1:ne,nd+1:end));
FF = rr(1:ncont,:)*Fsf;
XX = InfoEst.Xopt;
szXX = size(XX);
XX = reshape(XX,[szXX(1:2),size(Lest.Domain)]);
XX = pmat(XX,Lest.Domain);
trXXFF = trace( FF*(XX\FF'));
GamEst = sqrt( trXXFF);
GamEst = lpvmax(GamEst);

% Store Results
Gamma = sqrt(GamSF^2+GamEst^2);
Info.GAM = [GamSF GamEst];
Info.SF = InfoSF;
Info.Est = InfoEst;
Info.Fsf = Fsf;
Info.Lest = Lest;

