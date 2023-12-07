%% Simulink Blocks in LPVTools
%
% LPVTools provides Simulink blocks to interface to the state-space 
% LPV objects: |pss|, |upss| and |plftss|. 
% The Simulink blocks enable users to include LPV systems in Simulink
% simulation models.
% One Simlink block is for systems
% that depend on a time-varying parameter and its derivative, as seen in Equation 1, 
% while the other is for systems that do not depend explicitly on the
% derivative, as seen in Equation 2.
% 
% $$\left[ \begin{array}{c} \dot x (t) \\ y (t)\end{array} \right]
% = \left[ \begin{array}{cc} A(\rho(t),\dot{\rho}(t)) & 
% B(\rho(t),\dot{\rho}(t)) \\ C(\rho(t),\dot{\rho}(t)) & D(\rho(t),\dot{\rho}(t))
% \end{array} \right] \left[ \begin{array}{c} x (t) \\ u (t)\end{array} \right]
% \ \ \ \ \ \ \ (1)$$
% 
% $$\left[ \begin{array}{c} \dot x (t) \\ y (t)\end{array} \right]
% = \left[ \begin{array}{cc} A(\rho(t)) & B(\rho(t)) \\ C(\rho(t)) & D(\rho(t))
% \end{array} \right] \left[ \begin{array}{c} x (t) \\ u (t)\end{array} \right]
% \ \ \ \ \ \ \ (2)$$
% 
% The Simulink block for the system in Equation 2 is shown in Figure 1.
%
% <<LPVBlock.png>>
%
% _Figure 1: Simulink LPV block and block mask._
% 
% The block in Figure 1 has inputs for the system
% input $u(t)$ and the parameter vector $\rho(t)$, and an output for $y(t)$. 
% The block mask
% contains entries for the user to specify the system variable name, the
% order of the input parameter vectors $\rho$, and the state initial condition $x(0)$.
% The block is implemented as a Simulink S-function under the
% block mask.  The block currently performs a multidimensional linear
% interpolation to evaluate the state-space matrices at the specified
% parameter vector.  An efficient implementation of this linear
% interpolation has been coded to reduce computation and speed up the
% simulation time.
% 
% LPVTools also includes a block for systems that depend explicitly on both the
% time-varying parameter and its derivative, as seen in Equation 1. 
% This block is shown in Figure 2.
% 
% <<LPVBlock_withRates.png>>
% 
% _Figure 2: Rate-dependent Simulink LPV block and block mask._
