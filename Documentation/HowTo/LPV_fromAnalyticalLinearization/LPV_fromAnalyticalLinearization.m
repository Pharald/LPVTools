%% Deriving LPV models from Analytical Jacobian linearization
% 
%% Introduction
% 
% Consider a nonlinear system: 
% 
% $$ \dot x (t) = f(x(t),u(t),\rho(t)) \ \ \ \ \ \ \ (1)$$
% 
% $$ y(t) = h(x(t), u(t),\rho(t)) \ \ \ \ \ \ \ (2)$$
% 
% Where $x \in \mathcal{R}^{n_x}$, $y \in \mathcal{R}^{n_y}$, 
% $u \in \mathcal{R}^{n_u}$, and $\rho \in \mathcal{R}^{n_{rho}}$.
%
% Assume that $\rho(t) = \rho_0$ and $u(t) = \bar{u}(\rho_0)$ are 
% constant $\forall t \geq 0$.
% Then the solution of the nonlinear system is any 
% $x(t) = \bar{x}(\rho_0)$ and $y(t) = \bar{y}(\rho_0)$, such that 
% if $x(0) = \bar{x}(\rho_0)$, then $\forall t \geq 0$:
%
% $$ \dot{x}(t) = 0 = f(\bar{x}(\rho_0),\bar{u}(\rho_0),\rho_0)
% \ \ \ \ \ \ \ (3)$$
% 
% $$ y(t) = \bar{y}(\rho_0) = h(\bar{x}(\rho_0),\bar{u}(\rho_0),\rho_0)
% \ \ \ \ \ \ \ (4)$$
% 
% When $\rho(t)$ is a function of time, then the equilibrium 
% $\left( \bar{x}(\rho(t)), \bar{u}(\rho(t)), \bar{y}(\rho(t))\right)$
% is, in general, not a solution of the nonlinear system:
%
% $$ \frac{d}{dt}\bar{x}(\rho(t)) \neq 0 =
% f(\bar{x}(\rho(t)),\bar{u}(\rho(t)),\rho(t)) \ \ \ \ \ \ \ (5)$$
%
% We can linearize around 
% $\left( \bar{x}(\rho(t)), \bar{u}(\rho(t)), \bar{y}(\rho(t))\right)$
% even though it is not, in general, a solution of the nonlinear system. 
% Lets define perturbed quantities:
%
% $$\delta_x(t) = x(t) - \bar{x}(\rho(t))\ \ \ \ \ \ \ (6)$$
% 
% $$\delta_u(t) = u(t) - \bar{u}(\rho(t))\ \ \ \ \ \ \ (7)$$
% 
% $$\delta_y(t) = y(t) - \bar{y}(\rho(t))\ \ \ \ \ \ \ (8)$$
% 
% Using Taylor series expansion about 
% $\left( \bar{x}(\rho(t)), \bar{u}(\rho(t)), \bar{y}(\rho(t))\right)$,
% the system dynamics can be expressed as
% (dropping the notational dependence on time):
% 
% $$f(x,u,\rho) = f(\bar{x}(\rho),\bar{u}(\rho),\rho)+A(\rho)\delta_x
% +B(\rho)\delta_u + \Delta_f(\delta_x,\delta_u,\rho)
% \ \ \ \ \ \ \ (9)$$
% 
% $$h(x,u,\rho) = h(\bar{x}(\rho),\bar{u}(\rho),\rho)+C(\rho)\delta_x
% +D(\rho)\delta_u + \Delta_h(\delta_x,\delta_u,\rho)
% \ \ \ \ \ \ \ (10)$$
% 
% where $f(\bar{x}(\rho),\bar{u}(\rho),\rho) = 0$,
% $h(\bar{x}(\rho),\bar{u}(\rho),\rho) = \bar{y}(\rho)$,
% $\Delta_f(\delta_x,\delta_u,\rho)$ and
% $\Delta_h(\delta_x,\delta_u,\rho)$ terms represent higher-order terms
% of the Taylor series approximations, and
% 
% $$A(\rho) = \frac{\partial}{\partial x}f(x,u,\rho)
% \bigg|_{(x,u)=\left(\bar{x}(\rho),\bar{u}(\rho)\right)}
% \ \ \ \ \ \ \ (11)$$
%
% $$B(\rho) = \frac{\partial}{\partial u}f(x,u,\rho)
% \bigg|_{(x,u)=\left(\bar{x}(\rho),\bar{u}(\rho)\right)}
% \ \ \ \ \ \ \ (12)$$
%
% $$C(\rho) = \frac{\partial}{\partial x}h(x,u,\rho)
% \bigg|_{(x,u)=\left(\bar{x}(\rho),\bar{u}(\rho)\right)}
% \ \ \ \ \ \ \ (13)$$
%
% $$D(\rho) = \frac{\partial}{\partial u}h(x,u,\rho)
% \bigg|_{(x,u)=\left(\bar{x}(\rho),\bar{u}(\rho)\right)}
% \ \ \ \ \ \ \ (14)$$
% 
% Using this Taylor series approximation to linearize the dynamics of the 
% nonlinear system, yields:
% 
% $$ \begin{array}{l@{}l}
%  \frac{d}{dt}\delta_x  
% &{}= \frac{d}{dt}\left( x - \bar{x}(\rho)\right) \\
% &{}= \dot{x} - \frac{d}{dt}\bar{x}(\rho) \\
% &{}= f(x,u,\rho) - \frac{d}{dt}\bar{x}(\rho) \\
% &{}= A(\rho)\delta_x +B(\rho)\delta_x+
% \Delta_f(\delta_x,\delta_u,\rho)- \frac{d}{dt}\bar{x}(\rho) \end{array}
% \ \ \ \ \ \ \ (15)$$
% 
% 
%
% Similarly, the Taylor series approximation of $h(x,y,\rho)$ can be used
% to linearize the output equation:
%
% $$\begin{array}{l@{}l} \delta_y
% &{}= y-\bar{y}(\rho)\\
% &{}= h(x,y,\rho)-\bar{y}(\rho)\\
% &{}= \left[ \bar{y}(\rho)+C(\rho)\delta_x +D(\rho)\delta_x+
% \Delta_h(\delta_x,\delta_u,\rho)\right] -\bar{y}(\rho)\\
% &{}= C(\rho)\delta_x +D(\rho)\delta_x+\Delta_h(\delta_x,\delta_u,\rho)
% \end{array}
% \ \ \ \ \ \ \ (16)$$
%
% The final LPV model is thus:
% 
% $$ \begin{array}{l@{}l}
% \frac{d}{dt}\delta_x 
% &{}=A(\rho)\delta_x +B(\rho)\delta_x+
% \Delta_f(\delta_x,\delta_u,\rho)- \frac{d}{dt}\bar{x}(\rho) \\
% \delta_y &{}= C(\rho)\delta_x +D(\rho)\delta_x+\Delta_h(\delta_x,\delta_u,\rho)
% \end{array}
% \ \ \ \ \ \ \ (17)$$
%
% *Approximations*
% 
% Standard LPV approach is to neglect higher order terms $\Delta_f$ and
% $\Delta_h$, and the $-\dot{\bar{x}}$ term. However, the $-\dot{\bar{x}}$
% term can be retained and treated as a measurable disturbance. This
% can be expressed as $-\dot{\bar{x}} = G(\rho)\dot{\rho}$, 
% where $G(\rho) = -\frac{d\bar{x}(\rho)}{d\rho}$ 
% The higher order terms $\Delta_f$ and $\Delta_h$ are
% nonlinear functions. They can be handled (locally) as uncertainties using
% integral quadratic constraints. 
% 
% An interested reader, can refer to the work by Takarics and Seiler [1] 
% for additional details on this approach. If an analytical linearization
% is not possible, an LPV model can be constructed using numerical
% linearization directly from a nonlinear model (e.g. a Simulink model).
% Refer to section XXX for details.


%% Example
% 
% Consider the nonlinear system (from [2]): 
%
% $$\left[ \begin{array}{c} \dot x_1 (t) \\ \dot x_2 (t)\end{array} \right]
% = \left[ \begin{array}{cc} -1 & 0 \\ 1 & 0
% \end{array} \right] \left[ \begin{array}{c} x_1 (t) \\ x_2 (t)\end{array} \right]
% + \left[ \begin{array}{c} 1 \\ 0 \end{array} \right] u
% + \left[ \begin{array}{c} 0 \\ -x_2|x_2|-10\end{array} \right]
% \ \ \ \ \ \ \ (18)$$
%
% $$y(t) = x_2 \ \ \ \ \ \ \ (19)$$
% 
% Lets assume that we are given the control objective to make the output 
% $y(t)$ track a reference command $r(t)$. We will frame this as a LPV
% control problem, and derive a LPV model of this nonlinear model for this
% problem.
% 
% Let the desired operating point be scheduled y $\rho = r$. In this 
% formulation neither the dynamics ($f$ in Equation (1)), nor the output equation
% ($h$ in Equation (2)) directly depend on $\rho$.
%
% The equilibrium point, parameterized by $\rho$ is given by:
%
% $$ \bar{x}_1(\rho) = \rho|\rho|+10 \ \ \ \ \ \ \ (20)$$
% 
% $$ \bar{x}_2(\rho) = \rho \ \ \ \ \ \ \ (21)$$
%
% $$ \bar{u}(\rho) = \rho|\rho|+10 = \bar{x}_1(\rho) \ \ \ \ \ \ \ (22)$$
% 
% $$\bar{y}(\rho) = \rho \ \ \ \ \ \ \ (23)$$
% 
% Applying the approach described above, the nonlinear system in
% Equations (18)-(19) is linearized about the parameterized equilibrium point
% to obtain a LPV system:
%
% $$\dot{\delta_x} = 
% \left[ \begin{array}{cc} -1 & 0 \\ 1 & -2|\rho|) \end{array} \right] \delta_x 
% + \left[ \begin{array}{c} 1 \\ 0 \end{array} \right] u
% + \left[ \begin{array}{c} -2|\rho| \\ 1 ) \end{array} \right]
% \ \ \ \ \ \ \ (24)$$
%
% $$\delta_y = \left[ \begin{array}{cc} 0 & 1 \end{array} \right] \delta_x
% \ \ \ \ \ \ \ (25)$$
% 
% By formulating the control problem in the form of a LPV system which
% described the behaviour of the nonlinear system about a desired reference
% command, we have recast the problem into a regulation problem: $\delta_y
% = y-\bar{y}(\rho) = y - \rho$ and the control objective is to regulate
% $\delta_y(t) \to 0$ in the LPV model.
% 
% If we neglect the $-\dot{\bar{x}}$ term in the LPV system of Equations (24)-(25), 
% a grid-based LPV model of the system for $\rho \in [-5 0 10]$ can be
% constructed using the following commands:

% Define the parameter
p = pgrid('p',[-5 0 10]);

% Define the system matrices
A = [-1 0;1 -2*abs(p)];
B = [1;0];
C = [0 1];

% Define the grid-based LPV model
sys = ss(A,B,C,0)

%% 
% If we treat the $-\dot{\bar{x}}$ term as a exogenous disturbance to the 
% model then the grid-based LPV system can be modeled as:
Bd = [-2*abs(p);-1];
sys_dis = ss(A,[B Bd],C,0)

%%
% The $\dot{\rho}$ term in $-\dot{\bar{x}}$ is now an input to the model. 
% It is being treated as an exogenous disturbance, that is independent of 
% $\rho$. This assumption is, in general, conservative.



%% References
% 
% # B. Takarics and P. Seiler, "Gain Scheduling for Nonlinear Systems 
% via Integral Quadratic Constraints," _accepted to the American Control
% Conference_, 2015.
% # D. J. Leith and W. E. Leithead, "Counter-Example to a Common LPV
% Gain-Scheduling Design Approach," _UKACC International Control 
% Conference_, 2000. 





