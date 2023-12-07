%% Quasi-LPV Systems
% 
% A quasi-LPV system refers to case where the parameter trajectory is
% a function of the system state and inputs, i.e. $\rho(x,u,t)$
% where $x$ is the system state, $u$ is the system input, and $t$ is time.  
% The term "quasi" is used to indicate that the parameter dependence 
% actually introduces a nonlinearity in the system dynamics. 
%
% The current implementation of LPVtools treats quasi-LPV systems just like  
% any other LPV system for the purposes of analysis and synthesis, 
% i.e. the parameter is assumed to vary independently of 
% the system's state and input. This assumption will result in conservative 
% results when working with quasi-LPV systems in LPVTools. However,
% quasi-LPV systems can be simulated correctly using the |lpvlsim|,
% |lpvstep|, |lpvimpulse|, and |lpvinitial| functions in LPVTools.
% Specifically, the the parameter trajectory that is used in the simulation
% can be specified as a function of the state, input, and time.



