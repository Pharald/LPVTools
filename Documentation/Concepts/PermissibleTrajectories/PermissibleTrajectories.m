%% Permissible Parameter Trajectories


%%  Introduction
% 
% An LPV system is a time-varying, state-space model of the form:
%
% $$\left[ \begin{array}{c} \dot x (t) \\ y (t)\end{array} \right]
% = \left[ \begin{array}{cc} A(\rho(t)) & B(\rho(t)) \\ C(\rho(t)) & D(\rho(t))
% \end{array} \right] \left[ \begin{array}{c} x (t) \\ u (t)\end{array} \right]
% \ \ \ \ \ \ \ (1)$$
% 
% The LPV model in Equation (1) describes how the LPV system depends on a 
% set of time-varying parameters. 
% Its important to understand that for practical applications 
% (e.g. analysis in the LPV framework) each time-varying parameter in (1) has
% associated with it a set of _permissible parameter trajectories_, 
% which describe how the parameter can change with time in the model.
% The permissible parameter trajectories contstrain the parameter values to
% those for which the model is valid.
% 
% The set of allowable trajectories for a particular parameter satisfies 
% two properties: First, the parameter's value remains inside some interval
% of allowable values $[\rho_{min}, \rho_{max}]$ (an interval on the real line). 
% Second, the parameter's rate of change lies inside some interval
% $[\overline{\nu}, \underline{\nu}]$ (also an interval on the real line).  
% Hence, for an LPV system that only depends on 
% a single parameter $\rho \in \mathcal{R}$, a permissible trajectory is 
% any trajectory such that: $\rho_{min} \leq \rho(t) \leq \rho_{max}$
% and $\underline{\nu} \leq \dot{\rho}(t) \leq \overline{\nu}$ for all $t$.
% A trajectory is said to be "rate unbounded" if $\overline{\nu} = \infty$ 
% and $\underline{\nu} = -\infty$.
% 
% 
%% Example
% 
% Lets assume the LPV model in Equation (1) represents an aircraft, 
% and that the model is scheduled on altitude $h$. 
% If this particular model is only valid for a limited range of
% altitudes $h$: $5000~ft \leq h \leq 10000~ft$, and for slow
% variations in altitude $-10~ft/sec \leq \dot{h} \leq 10~ft/sec$, then 
% the set of permissible parameter trajectories contains any altitude
% trajectory such that 
% 
% $$h(t) \in [5000,10000]~ft,~\forall t$$  
%  
% and  
% 
% $$\dot{h}(t) \in [-10,10]~ft/sec,~\forall t$. 
%
% An example of a permissible parameter trajectory for this system 
% is shown in Figure 1.
%
% <<AltitudeExample.png>>
%
% _Figure 1: A permissible altitude trajectory._
%
% 
%% Formal Definition 
% 
% Given an LPV system that depend on a set of time-varying parameters 
% $\rho \in \mathcal{R}^{n_\rho}$. A permissible parameter 
% trajectory is any trajectory such that $\rho$ lies inside
% the compact set $\mathcal{P} \subseteq \mathcal{R}^{n_\rho}$
% and  $\dot{\rho}(t)$ lies inside the set 
% $\mathcal{D} \subseteq \mathcal{R}^{n_\rho}$. 
% The set $\mathcal{P}$ is the $n_\rho$ dimensional hyper rectangle formed  
% by $[\rho_{1,min},\rho_{1,max}]\times[\rho_{2,min},\rho_{2,max}]\times
% \ldots \times [\rho_{n_\rho,min},\rho_{n_\rho,max}]$, 
% and the set $\mathcal{D}$ is the $n_\rho$ dimensional hyper rectangle formed 
% by $[\underline{\nu}_1,\overline{\nu}_1]\times[\underline{\nu}_2,\overline{\nu}_2]\times
% \ldots \times [\underline{\nu}_{n_\rho},\overline{\nu}_{n_\rho}]$. 



