%% Simulation of LPV systems
% A primer on simulation in LPVTools.

%% Introduction
% 
% LPVTools provides a set of tools to perform time-domain simulation of LPV
% systems. The tools can be split into two parts: First, there are 
% overloaded versions of Linear Time-Invariant (LTI) simulation tools from 
% the Control Systems Toolbox (|lsim|, |step|, |impulse|, |initial|). 
% These functions can be used to evaluate the pointwise behaviour of the 
% LPV system when the scheduling parameter ($\rho$) is held constant. 
% Second, there are LPV simulation tools that enable simulation of the
% LPV system when the parameter is allowed to vary with time. These LPV
% simulation tools are able to capture the time-varying nature of the 
% LPV system's dynamics, and allow the system's behaviour to be evaluated
% for different parameter trajectories.

%% LPV Simulation Commands
% 
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVLSIM\html\LPVLSIMdoc.html">LPVLSIM</a> </td>
% <td>Simulate parameter dependent time response of a LPV system.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSTEP\html\LPVSTEPdoc.html">LPVSTEP</a></td>
% <td>Simulate parameter dependent step response of a LPV system.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVINITIAL\html\LPVINITIALdoc.html">LPVINITIAL</a> </td>
% <td>Simulate initial conditions response of a LPV system.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVIMPULSE\html\LPVIMPULSEdoc.html">LPVIMPULSE</a></td>
% <td>Simulate parameter dependent impulse response of a LPV system.</td>
% </tr>
% </table>
% </html>


%% Examples and How To
% 
% * <..\..\HowTo\SimulateQuasiLPV\html\SimulateQuasiLPV.html Tutorial: Simulate a Quasi-LPV System>
% * <..\..\GettingStarted\MiniTutorials\Grid_Analysis\html\Grid_Analysis.html Tutorial: Analysis of a grid-based LPV system>
% * <..\..\GettingStarted\MiniTutorials\Grid_Synthesis\html\Grid_Synthesis.html Tutorial: Control of a grid-based LPV system>
% * <..\..\Demos\SpinningDisk_Stochastic\html\SpinningDisk_Stochastic.html Example: Stochastic LPV control of spinning mass>
% * <..\..\Demos\SpinningDisk_L2_LFT\html\SpinningDisk_L2_LFT.html Example: LPV control of spinning mass using LFT framework>


%% Concepts
% 
% * <..\..\Concepts\PermissibleTrajectories\html\PermissibleTrajectories.html Permissible
% Parameter Trajectories>
% * <..\..\Concepts\QuasiLPV\html\QuasiLPV.html Quasi-LPV Models>
