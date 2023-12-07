%% Control Design for Stochastic LPV Systems
% 
%  The following example from Wu [1] describes a control design for a
%  mechanical system in the stochastic LPV framework. 
%
%% Problem Statement: Coupled Spinning Disks
%
% A pair of rotating disks is shown in Figure 1. The slotted disks
% rotate at rates $\Omega_1$ rad/s and $\Omega_2$ rad/s. The disks contain
% masses $M_1$ and $M_2$, which can move in the horizontal plane, and slide 
% in the radial direction The two masses are connected by a wire, which can
% transmitt force in both compression and tension. The wire acts as a
% spring with spring constant $k$. A coupler and pulley system handles the
% wire routing between the two disks. The effect of friction due to the 
% sliding motion of each mass is modeled by a damping 
% coefficient $b$.
%
% <<SpinningDisks.png>>
%
% _Figure 1: A pair of spinning disk with masses $M_1$ and $M_2$ [1]._
% 
% The problem is to control the position of the two masses: $r_1$ for mass
% $M_1$ and $r_2$ for $M_2$. The control input is a radial force $f_u$
% acting on mass $M_1$, while there is a radial disturbance force
% acting on each mass: $f_1$ acting on $M_1$ and $f_2$ acting on $M_2$.
% The equations of motion for this problem are:
% 
% $$
% \begin{array}{l@{}l}
% M_1 \left(\ddot{r}_1-\Omega_1^2 r_1 \right) + b \dot{r}_1 + k(r_1 +r_2) 
% &{}= f_u + f_1 \\
% M_2 \left(\ddot{r}_2-\Omega_2^2 r_2 \right) + b \dot{r}_2 + k(r_1 +r_2) 
% &{}= f_2 \end{array}
% \ \ \ \ \ \ \ (1)$$
%
% where: 
% 
% * $M_1 = 1$ kg is the mass of the sliding mass in disk 1.
% * $M_2 = 0.5$ kg is the mass of the sliding mass in disk 2.
% * $b = 1$ kg is the damping coefficient due to fricition (kg/sec).
% * $k = 200$ N/m is the spring constant of the wire connecting the masses.
% * $r_1(t)$ is the position of the $M_1$ sliding mass relative to the center of disk 1 (m). 
% * $r_2(t)$ is the position of the $M_2$ sliding mass relative to the center of disk 2 (m). 
% * $\Omega_1(t)$ is the rotational rate of disk 1 (rad/s). 
% * $\Omega_2(t)$ is the rotational rate of disk 2 (rad/s). 
% * $f_u$ is the control force acting radially on mass $M_1$ (N).
% * $f_1$ is the disturbance force acting radially on mass $M_1$ (N).
% * $f_2$ is the disturbance force acting radially on mass $M_2$ (N).
% 
%
% The rotational rates of the spinning disk are allowed to vary: 
% $\Omega_1 \in \left[0, 3\right]$ rad/s and $\Omega_2 \in \left[0, 5\right]$
% The rotational rates are not known in advance but are measured and
% available for control design.
% 
% The objective of the control design is to command the radial 
% position of mass $M_2$. Note that the control input is applied to mass
% $M_1$, and is transmitted to mass $M_2$ through the wire connecting the 
% two disks. 


%% LPV System Formulation:
% 
% The system in Equation (1) is already processing the inputs and outputs
% linearly, and the only nonlinear elements in Equation (1) are the 
% rotation rates $\Omega_1$ and $\Omega_2$. Lets choose the rotation rates 
% as the scheduling parameters and transform the system into a LPV model.
% Define $\rho_1 = \Omega_1^2$ and $\rho_2 = \Omega_2^2$, such that 
% $\rho_1 \in [0,25]$ and $\rho_2 \in [0 9]$. Furthermore, let
%
% $$ 
% \left[ \begin{array}{c}  x_1 (t) \\ x_2 (t) \\
% x_3 (t) \\ x_4 (t) \end{array} \right] 
%  = 
% \left[ \begin{array}{c}  r_1 (t) \\ r_2 (t) \\
% \dot r_1 (t) \\ \dot r_2 (t) \end{array} \right]
% \ \ \ \ \ \ \ (2)$$
% 
% and 
% 
% $$
% \left[ \begin{array}{c}  u (t) \\ d_1 (t) \\ d_2 (t) \end{array} \right] 
%  = 
% \left[ \begin{array}{c} f_u (t) \\ f_1 (t) \\ f_2 (t) \end{array} \right]
% \ \ \ \ \ \ \ (3)$$
% 
% And rewrite the system in Equation (1) as the LPV system:
%
% $$
% \left[ \begin{array}{c} \dot x_1 (t) \\ \dot x_2 (t) \\
% \dot x_3 (t) \\ \dot x_4 (t) \end{array} \right]
%  = \left[ \begin{array}{cccc} 
% 0 & 0 & 1 & 0 \\ 
% 0 & 0 & 0 & 1 \\ 
% \rho_1 - \frac{k}{M_1} & -\frac{k}{M_1} & -\frac{b}{M_1} & 0 \\
% -\frac{k}{M_2} & \rho_2 - \frac{k}{M_2} & 0 & -\frac{b}{M_2} \\ 
% \end{array} \right]
% \left[ \begin{array}{c}  x_1 (t) \\ x_2 (t) \\
% x_3 (t) \\ x_4 (t) \end{array} \right]
% +
% \left[ \begin{array}{ccc} 
% 0 & 0 & 0 \\ 
% 0 & 0 & 0  \\ 
% \frac{1}{M_1} & \frac{0.1}{M_1} & 0 \\
% 0 & 0 & \frac{0.1}{M_2} \\ 
% \end{array} \right]
% \left[ \begin{array}{c}  f_u (t) \\ f_1 (t) \\ f_2 (t)\end{array} \right]
% \ \ \ \ \ \ \ (4)$$
%
% $$ 
% y(t) = \left[\begin{array}{cccc} 0 & 1 & 0 & 0 \end{array} \right]
% \left[ \begin{array}{c}  x_1 (t) \\ x_2 (t) \\
% x_3 (t) \\ x_4 (t) \end{array} \right]
% \ \ \ \ \ \ \ (5)$$
% 
% The following commands will generate the LPV system in Equations (5)-(6) 

% Define system parameters
m1 = 1;
m2 = 0.5;
k = 200;
b = 1;

% Define timve-varying parameters
rho1 = pgrid('rho1',[0 4.5 9]);
rho2 = pgrid('rho2',[0 12.5 25]);

% Define system matrices:
A = [   0              0            1       0    ; ...
        0              0            0       1    ; ...
    rho1-k/m1    -k/m1        -b/m1     0    ; ...
      -k/m2         rho2-k/m2   0     -b/m2  ];

B = [    0      0       0     ; ...
         0      0       0     ; ...
      1/m1   .1/m1      0     ; ...  
         0      0     .1/m2 ];

C = [ 0 1 0 0];

D = [0 0 0];

% Define the LPV system:
sys = ss(A,B,C,D)


%% LPV Control Design
%
% We will design an LPV controller that optimizes the stochastic bound 
% 
% <<IC.png>>
%
% _Figure 2: The weighted design interconnection [1]._
% 
% The weighted interconnection shown in Figure 2 will be used for
% synthesis. With appropriate weights selected, this interconnection, 
% defines the desired performance of the LPV control law. The controller
% will be synthesized to yield optimized performance in terms of this
% specification (additional details on formulating design problems in this 
% way are provided in the Robust Control Toolbox's documentation on
% <matlab:open(fullfile(docroot,'robust/gs/h-infinity-performance.html')) 
% characterizing closed-loop performance objectives>.)
% 
% 
% The weights for the interconnection are chosen as:

% Define weights:
We = tf(2,[1 0.04]);
Wu = 1/50;
Wa = 0.00001;
Wn = tf([1,0.4],[0.01 400]);
act = ss(-100,100,1,0);

%%
% and we can form the weighted interconnection shown in Figure 2:

% Form synthesis interconnection:
systemnames = 'sys act We Wu Wa Wn';
inputvar = '[ref;d{2};dn;u]';
outputvar = '[Wu;Wa;We;ref;sys+Wn]';
input_to_sys = '[act;d]';
input_to_act = '[u]';
input_to_Wa = '[act]';
input_to_Wu = '[u]';
input_to_Wn = '[dn]';
input_to_We = '[sys-ref]';
G = sysic;

%%
% We will use the function |lpvstochsyn| to synthesize a LPV controller
% which minimizes the 
% <..\..\Concepts\StabilityAndInducedGain\html\StabilityAndInducedGain.html stochastic LPV bound >
% of the weighted interconnection |G| in Figure 2. This synthesis will
% assume that the rate of variation in $\rho_1$ and $\rho_2$ are unbounded,
% i.e. there is no limit to fast the parameters can change with time.

% Synthesize two-degree of freedom controller.
nmeas = 2;					% # of measurements
ncont = 1;					% # of controls
[Knr,Gamma,Info] = lpvstochsyn(G,nmeas,ncont);

%%
% The control design is successfull, and the controller |Knr| is guarenteed 
% to achieve a stochastic LPV bound that is less than or equal to |Gamma|
Gamma

%% Computing a Tighter Bound on the Stochastic LPV Bound
% 
% The value reported in |Gamma| is an only an upper bound on the stochastic 
% LPV bound. 
% We will now compute a tighter bound on the stochastic LPV bound achieved  
% by |Knr|. We will insert |Knr| into the weighted interconnection and use
% this closed-loop weighted interconnection to
% compute the H2 norm for the Linear 
% Time-Invariant (LTI) system obtained by holding the parameters constant at each
% grid point in the domain of the LPV system: $\left[\rho_1,\rho_2 \right]
% = \left[ 0, 4.5, 9\right] \times\left[0, 12.5,25\right]$.
% The stochastic LPV bound will always be greater or equal to the largest 
% H2 norm on the grid of parameter values. Hence, the largest H2 norm we 
% compute will provide a lower bound on the stochastic LPV bound.

% Insert Knr into the weighted interconnection:
CLICnr = lft(G,Knr);

% Compute the H2 norm for the 
LTI_CLnorm = norm(CLICnr,2);

% extract the largest H2 norm on the grid of parameter values:
lpvmax(LTI_CLnorm)

%%
% The results indicate that the stochastic LPV bound achieved by |Knr| is
% between 1.5506 and 1.5998. 

%% Evaluating Pointwise Performance 
% 
% Lets look at the closed-loop response for the original LPV system in the
% loop with the controller |Knr|:

% Form closed-loop sytem
systemnames = 'sys act Knr';
inputvar = '[r;d{2}]';
outputvar = '[sys]';
input_to_sys = '[act;d]';
input_to_act = '[Knr]';
input_to_Knr = '[r;sys]';
CL = sysic;

%%
% We will start by applying LTI analysis techniques to evaluate the
% performance of the LPV controller. We will evalaluate the LTI closed-loop
% system obtained at each grid point in the domain $\left[\rho_1,\rho_2 \right]
% = \left[ 0, 4.5, 9\right] \times\left[0, 12.5,25\right]$, when the
% parameter is held constant.
% 
% We can start by studying the pointwise step response. The Control System
% Toolbox's command |step| is overloaded to work with the |pss| object. It
% will compute the step response at each point in the grid, and plot them
% together on the same plot:

step(CL)

%%
% The controller achives good tracking with minimal overshoot, 
% approximatly 2.5 sec settling time and less than 3% steady-state tracking
% error. Furthermore, the effect of the disturbances is minimal.
% 
% Similarly we can create pointwise frequency responses using the |bode| 
% command. Lets compare the pointwise frequency response of the open- and
% closed-loop system at frequencies between 0.1 to 20 rad/s

bode(sys,{0.1,100},'b')
hold on
bode(CL,{0.1,100},'r--')
legend('Open-loop','Closed-loop','location','best')
grid minor
hold off

%%
% The closed-loop system has a bandwidth of approximatly 3 rad/s.

%% Evaluating LPV Time-Domain Performance 
% 
% We will use the LPV simulation capabilities to inspect the performance
% of |Knr|. The parameters $\rho_1$ and $\rho_2$ will be made to vary with 
% time: $\rho_1(t) = sin(t)+1.5$ and $\rho_2(t) = 0.5*cos(5*t)+3$, while
% the system tracks a unit step response and rejects disturbances 
% $d_1 = cos(3t)+sin(5t)+cos(8t)$ and $d_2 = cos(t)+sin(2t)+cos(4t)$.

% Define the inputs to the system:
t = [0:0.01:15]';
u = ones(size(t));
d1 = cos(3*t)+sin(5*t)+cos(8*t);
d2 = cos(t)+sin(2*t)+cos(4*t);


% Define the trajectories of the parameters:
ptraj.time = t;
ptraj.rho1 = sin(t)+1.5;
ptraj.rho2 = .5*cos(5*t)+3;

% Perform LPV simulation:
lpvlsim(CL,ptraj,[u,d1,d2],t);


%%
% Lets evaluate the impact that the disturbances had on the response by
% repeating the simulation without disturbances:

lpvlsim(CL,ptraj,[u,d1,d2],t);
hold on
lpvstep(CL(:,1),ptraj);
legend('LPV simulation with disturbances',...
       'LPV simulation without disturbances',...
       'location','best')



%%
% The error due to the 
% distrubances is on the order of 5%, and they are confined to an
% approximately 1-2 Hz oscillation about the nominal. This frequency range
% is incidentally where the system's pointwise frequency response analysis 
% indicated that disturbances would have the greatest effect on the output 
% signal.


%% References
%
% # F. Wu, "Control of Linear Parameter Varying Systems," Ph.D. dissertation, 
% Department of Mechanical Engineering, University of California at 
% Berkeley, CA, May 1995.

