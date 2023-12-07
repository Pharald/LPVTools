%% LPVMIXSYN - Parameter-varying mixed-sensitivity synthesis
%
%  
%% Syntax
%
%    [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3)
%    [K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,Xb,Yb)
%    [K,GAM,INFO]=lpvmixsyn(G,...,OPT)
%
%% Description
% 
% |[K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3)| synthesizes the parameter-varying
% mixed-sensitivity controller |K|, which minimizes the induced $L_2$ norm of  
% |W1*S|, |W2*K*S| and |W3*T|, where |S = inv(I+G*K)|, |T = G*K*inv(I+G*K)|, and 
% |W1|, |W2| and |W3| are stable |pss|, |pmat| or |double| weights of appropriate size. 
% |GAM| is the induced $L_2$ norm acheived by |K|. |INFO| is a structure containing 
% data from the Linear Matrix Inequalities that are solved to obtain |K|. 
% A call to |lpvmixsyn| without a |basis| function argument generates a 
% controller assuming no bounds on the parameter rate of variation. 
%
% |[K,GAM,INFO]=lpvmixsyn(G,W1,W2,W3,Xb,Yb)| computes the rate-bounded 
% mixed-sensitivity controller |K|, where the rate-bounds of the independent 
% variables of |G|, |W1|, |W2| and |W3| are incuded in the synthesis conditions. 
% |Xb| and |Yb| are |basis| objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for |K|.
%
% |[K,GAM,INFO]=lpvmixsyn(G,...,OPT)| allows the user to pass in a 
% |lpvsynOptions| object. 
