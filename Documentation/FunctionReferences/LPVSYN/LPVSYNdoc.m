%% LPVSYN - Parameter-dependent controller synthesis for LPV systems
%
%  
%% Syntax
%
%    [K,GAM,INFO] = lpvsyn(P,NMEAS,NCON)
%    [K,GAM,INFO] = lpvsyn(P,NMEAS,NCON,Xb,Yb)
%    [K,GAM,INFO] = lpvsyn(P,NMEAS,NCON,...,OPT)
% 
%% Description
% 
% |[K,GAM,INFO] = lpvsyn(P,NMEAS,NCON)| computes a parameter-varying
% controller |K| which minimizes the induced $L_2$ norm of the interconnection 
% defined by |lft(P,K)|. |K| is a |pss| or |plftss| with |NMEAS| inputs and |NCON| outputs, 
% defined on same domain as |P|. |GAM| is the induced $L_2$ norm of |lft(P,K)|.
% This three argument call assumes that the rate-bounds of the independent
% variables in |P| are |[-inf,inf]|. |INFO| is a structure containing data from
% the Linear Matrix Inequalities that are solved to obtain |K|.
%
% |[K,GAM,INFO] = lpvsyn(P,NMEAS,NCON,Xb,Yb)| computes the rate-bounded 
% parameter-varying controller |K| for a system |P|. |K| is the controller which 
% minimizes the induced $L_2$ norm of |lft(P,K)| when the rate-bounds of the  
% independent variables of |P| are incorporated into the synthesis. 
% |Xb| and |Yb| are |basis| objects, which describe the assumed parameter 
% dependence of the lyapunov matrices used in solving for |K|.
%
% |[K,GAM,INFO] = lpvsyn(P,NMEAS,NCON,...,OPT)| allows the user to pass in
% a |lpvsynoptions| object. 
%
% The default algorithm for |lpvsyn| will solve the given synthesis problem
% twice. The first iteration attempts to find a solution that minimizes the
% induced $L_2$ norm of |lft(P,K)|. The second iteration will solve the 
% optimization problem again, with the caveat that any solution that is 
% % found to lie within 15% of the optimal induced $L_2$ norm of |lft(P,K)| from 
% the first iteration, is satisfactory. This formulation has been found to 
% yield controllers that are better numerically conditioned. The back-off 
% factor of 15% can be reset to a different value in |lpvsynoptions|.
