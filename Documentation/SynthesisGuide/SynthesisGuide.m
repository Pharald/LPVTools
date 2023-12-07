%% Synthesis for LPV systems
% A primer on synthesis in the LPV framework.

%% LPV Synthesis
%
% *Problem Statement*
% 
% Consider a parameter dependent linear plant $P_\rho$ of the form:
% 
% $$\left[ \begin{array}{c} \dot x (t) \\ e (t) \\ y (t) \end{array} \right]
% = \left[ \begin{array}{ccc} A(\rho(t)) & B_1(\rho(t)) & B_2(\rho(t)) \\ 
% C_1(\rho(t)) & D_{11}(\rho(t)) & D_{12}(\rho(t)) \\
% C_2(\rho(t)) & D_{21}(\rho(t)) & D_{22}(\rho(t)) \end{array} \right] 
% \left[ \begin{array}{c} x (t) \\ d (t) \\ u (t) \end{array} \right]
% \ \ \ \ \ \ \ (1)$$
%
% where $\rho$ is a time varying parameter, that takes on values in a known
% compact set $\mathcal{P}$ and has known bound on $\dot{\rho}$, 
% $\overline{\nu} \leq \dot{\rho} \leq \underline{\nu}$. 
% The time variations of $\rho(t)$ are not known in advance, but 
% the parameter values are measured in real-time and available for
% control design.
%
% <<ClosedLoop.png>>
%
% _Figure 1: Closed-loop interconnection for LPV synthesis problem._
% 
% The control problem is to synthesize a controller $K_{\rho}$ 
% such that the closed-loop system shown in Figure 1, is stable and
% the gain from $d$ to $e$ is minimized.
% This requires that the controller be designed such that the closed-loop 
% performance is optimized in the presence of rate-bounded,
% time-varying parameter trajectories $\rho \in \mathcal{P} \subset \mathcal{R}^n_{\rho}$.
% Denote the closed-loop system by $lft(P_{\rho},K_{\rho})$, and the gain of
% this closed-loop system by $\|lft(P_{\rho},K_{\rho})\|$
% Then the design objective can be stated as:
%
% $$ \min_{K_{\rho}} \max_{\rho \in \mathcal{P}, 
%    \overline{\nu} \leq \dot{\rho} \leq \underline{\nu}} 
%    \|lft(P_{\rho},K_{\rho})\|
% \ \ \ \ \ \ \ (2)$$
% 
% The resulting controller is itself parameter dependent - using the
% available real-time information of the parameter variation. 
% In the grid-based LPV framework 
% 
% *LPVTools Implementation*
%
% LPVTools implements LPV controller synthesis for both the LFT-based
% LPV framework and the grid-based LPV framework.  
% The synthesis functions generate controllers which 
% optimize the performance of the
% closed-loop system while taking into account the 
% permissible parameter trajectories: $\rho \in \mathcal{P}$, 
% subject to $\overline{\nu} \leq \dot{\rho} \leq \underline{\nu}$.
%
% In the grid-based LPV framework |lpvsyn|, |lpvncfyn|, |lpvmixsyn|, 
% |lpvloopshape|, and |lpvstochsyn| are used to synthesize LPV 
% output-feedback controllers. 
% |lpvsfsyn| is used to synthesize LPV state-feedback controllers, 
% and |lpvestsyn| is used to generate LPV estimators.
% These functions can be used to generate controllers and estimators to 
% minimize either the induced $L_2$ norm  (based on results by Becker [1] and Wu [2,3], 
% with pole-constrained synthesis based on the derivation by Lee [4]) or the 
% stochastic LPV bound (based on results by Wu [2]).
% In the LFT-based LPV framework only |lpvsyn| is provided to 
% synthesize LPV output-feedback controllers, and it 
% implements an algorithm which minimizes the 
% induced $L_2$ norm (based on results by Packard [5], and 
% Apkarian and Gahinet [6]).
%
% The LPV controller synthesis conditions lead to a set of Linear Matrix
% Inequalities (LMIs) which must be solved in order to generate a controller.
% These LMIs suffer from similar computational issues to the
% <..\..\Concepts\StabilityAndInducedGain\html\StabilityAndInducedGain.html
% LPV analysis conditions>, and their complexity also grows with
% $O(2^{n_\rho})$. 
% 
% *References*
%
% # G. Becker, "Quadratic Stability and Performance of Linear Parameter 
% Dependent Systems," Ph.D. Dissertation, University of California,
% Berkeley, 1993.
% # F. Wu, "Control of Linear Parameter Varying Systems," PhD thesis, University of California,
% Berkeley, 1993.
% # F. Wu, X. Yang, A. Packard, and G. Becker, "Induced L2 norm control for LPV systems with
% bounded parameter variation rates," _International Journal of Nonlinear and Robust Control_,
% vol. 6, pp. 983-998, 1996.
% # L. H. Lee, "Identification and Robust Control of Linear 
% Parameter-Varying Systems," Ph.D. Dissertation, University of California 
% at Berkeley, 1997, doi:10.1.1.55.2269.
% # A. Packard, "Gain scheduling via linear fractional transformations," _Systems and Control
% Letters_, vol. 22, no. 2, pp. 79-92, 1994.
% # P. Apkarian and P.Gahinet, "A convex characterization of gain-scheduled Hinfinity controllers,"
% _IEEE Transactions on Automatic Control_, vol. 40, no. 5, pp. 853-864, 1995.
					 	                                  

%% LPV Synthesis Commands
% 
% LPVTools provides the following functions to design controllers for multiinput-multioutput (MIMO) LPV models:
% 
% 
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSYN\html\LPVSYNdoc.html">LPVSYN</a> </td>
% <td>Synthesize a LPV controller</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVNCFSYN\html\LPVNCFSYNdoc.html">LPVNCFSYN</a></td>
% <td>Normalized coprime factor LPV controller synthesis</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVLOOPSHAPE\html\LPVLOOPSHAPEdoc.html">LPVLOOPSHAPE</a></td>
% <td>LPV loop-shaping synthesis</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVMIXSYN\html\LPVMIXSYNdoc.html">LPVMIXSYN</a></td>
% <td>LPV mixed-sensitivity synthesis</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSFSYN\html\LPVSFSYNdoc.html">LPVSFSYN</a></td>
% <td>Synthesize a LPV state-feedback controller</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVESTSYN\html\LPVESTSYNdoc.html">LPVESTSYN</a></td>
% <td>Synthesize a LPV state estimator</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSTOCHSYN\html\LPVSTOCHSYNdoc.html">LPVSTOCHSYN</a></td>
% <td>Synthesize stochastic LPV controller</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSYNOPTIONS\html\LPVSYNOPTIONSdoc.html">LPVSYNOPTIONS</a></td>
% <td>Create options object for LPV synthesis and analysis</td>
% </tr>
% </table>
% </html>

%% Examples and How To
% 
% * <..\..\HowTo\BASISexample\html\BASISexample.html Tutorial: Creating basis functions>
% * <..\..\GettingStarted\MiniTutorials\Grid_Synthesis\html\Grid_Synthesis.html Tutorial: Synthesis for gridded LPV systems>
% * <..\..\GettingStarted\MiniTutorials\LFT_Tutorial\html\LFT_Tutorial.html Tutorial: Synthesis for LFT LPV systems>
% * <..\..\Demos\SpinningDisk_Stochastic\html\SpinningDisk_Stochastic.html Example: Stochastic LPV control of spinning mass>
% * <..\..\Demos\SpinningDisk_L2_LFT\html\SpinningDisk_L2_LFT.html Example: LPV control of spinning mass using LFT framework>

%% Concepts
% 
% * <..\..\Concepts\PermissibleTrajectories\html\PermissibleTrajectories.html Permissible
% Parameter Trajectories>
% * <..\..\Concepts\StabilityAndInducedGain\html\StabilityAndInducedGain.html Stability and Induced Gain>
% * <matlab:open(fullfile(docroot,'robust/gs/h-infinity-performance.html')) 
% Characterizing Closed-loop Performance Objectives>

%% LTI synthesis capabilities
% 
% Overloaded LTI synthesis function from the Control Systems Toolbox 
% and the Robust Control Toolbox are provided for LPV systems. 
% (e.g. |lqr|, |hinfsyn|, |h2syn|, |loopsyn|, |ncfsyn|, and |mixsyn|). 
% These functions perform the controller synthesis pointwise
% in the parameter domain of the controller.
% In each case the resulting controller is not a LPV controller 
% (i.e. one that satisfies the LPV analysis conditions), 
% but a collection of LTI controllers - one for each point.