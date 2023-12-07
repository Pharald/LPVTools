%% LPVNORM - Compute bound on norm for |pss| systems.
%
%  
%% Syntax
%
%    [Gamma,X] = lpvnorm(P)
%    [Gamma,X] = lpvnorm(P,'L2')
%    [Gamma,X] = lpvnorm(P,'LQG') 
%    [Gamma,X] = lpvnorm(P,Xb)
%    [Gamma,X] = lpvnorm(P,Xb,'L2')
%    [Gamma,X] = lpvnorm(P,Xb,'LQG')
%
%% Description
% 
% |lpvnorm| computes a bound on the norm of a |pss| system over the set of
% all permissible trajectories of the independend variables which the |pss| 
% depends on. 
%
% |[Gamma,X] = lpvnorm(P,'L2')| computes an upper bound |Gamma| on the induced 
% $L_2$ norm of the |pss| |P|. The upper bound |Gamma| and a constant (parameter
% independent) matrix |X| are computed to satisfy the induced norm linear
% matrix inequality (LMI) condition. |X| is returned as a |pmat|. The $L_2$ norm
% bound is valid for arbitrarily fast variations of the system parameters.
%
% |[Gamma,X] = lpvnorm(P,'LQG')| computes an upper bound |Gamma| on the 
% stochastic LPV bound. The stochastic LPV bound is defined as the expected 
% value of the average instantaneous power of the output of |P|, assuming its
% inputs are zero mean, white-noise processes with unit intensity.
%
% |[Gamma,X] = lpvnorm(P,Xb,ALG)| computes a tighter (less conservative) 
% bound on the norm by using a parameter dependent matrix |X| = $X(\rho)$ 
% and bounds on the parameter rates of variation. The basis functions used 
% to construct $X(\rho)$ are specified with the |basis| object |Xb|. |ALG| can be 
% either |'L2'| or |'LQG'|. A call without the |ALG| argument is equivalent to 
% |[Gamma,X] = lpvnorm(P,Xb,'L2')|.