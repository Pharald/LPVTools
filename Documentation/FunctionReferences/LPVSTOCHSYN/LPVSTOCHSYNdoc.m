%% LPVSTOCHSYN - LPV controller synthesis for stochastic LPV systems
%
%  
%% Syntax
%
%    [K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON)
%    [K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON,XB,YB)
%    [K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON,XB,YB,OPT)
% 
%% Description
%
% |[K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON)| computes a parameter-varying 
% controller for the |pss| |P|. |NCON| specifies 
% the number of available control inputs in |P|. |NMEAS| specifies the number 
% of available measurements being output from |P|. |K| is the parameter-dependent 
% controller for the stochastic plant |P|, which minimizes the stochastic LPV 
% bound |GAMMA|. The stochastic LPV bound is defined as the expected value of the 
% average instantaneous power of the output of |P|, assuming its inputs are 
% zero mean, white-noise processes with unit intensity. |INFO| is a struct 
% with additional data. 
% 
% |[K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON,XB,YB)| performs a rate-bounded 
% synthesis. |Xb| and |Yb| are |basis| objects which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for |K|.
%
% |[K,GAMMA,INFO] = lpvstochsyn(P,NMEAS,NCON,XB,YB,OPT)| allows the user to 
% pass in a |lpvsynoptions| object. 
%
% The default algorithm for |lpvstochsyn| will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% stochastic LPV bound of |lft(P,K)|. The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% found to lie within 15% of the optimal stochastic LPV bound of |lft(P,K)| 
% from the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be reset to a different value in |lpvsynoptions|.
