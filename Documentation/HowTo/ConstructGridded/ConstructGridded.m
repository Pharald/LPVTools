%% Constructing Gridded LPV Systems
% 
% 
%% LPV Models
%
% Linear Parameter-Varying Models
%
% *Functions*
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\..\FunctionReferences\LPVSYN\html\LPVSYNdoc.html">LPVSYN</a> </td>
% <td>Synthesize an LPV controller</td>
% </tr>
% <tr>
% <td><a href="..\..\..\FunctionReferences\LPVNCFSYN\html\LPVNCFSYNdoc.html">LPVNCFSYN</a></td>
% <td>Normalized coprime factor LPV controller synthesis</td>
% </tr>
% <tr>
% <td><a href="..\..\..\FunctionReferences\LPVSFSYN\html\LPVSFSYNdoc.html">LPVSFSYN</a></td>
% <td>Synthesize a LPV state-feedback controller</td>
% </tr>
% <tr>
% <td><a href="..\..\..\FunctionReferences\LPVESTSYN\html\LPVESTSYNdoc.html">LPVESTSYN</a></td>
% <td>Synthesize a LPV state estimator</td>
% </tr>
% </table>
% </html>

% 
%% Examples and How To

%% Concepts

Let $S(\rho)$ represent a state-space system of the form:

<LPVsystem.png>

which depends on a time-varying parameter vector $\rho \in \mathcal{R}^{n_\rho}$.
A grid-based LPV model of this system is a collection of linearizations 
on a gridded domain of parameter values.


The linearizations are obtained through Jacobian linearization at each 
grid point. Each linearization approximates the system's dynamics in the 
vicinity of a particular grid point, and the grid of linearizations 
captures the system's parameter dependence implicitly.
Hence, linearization based LPV models do not require any special 
dependence on the parameter vector. 

Jacobian Linearization is the predominant method of constructing grid-based LPV models.
In this case the nonlinear system to be approximated by a LPV model is linearized 
at each point in the desired grid of parameter values (e.g. batch linearization of Simulink models.)


Let $G$ be a nonlinear system that depends on Mach and altitude values. 



Figure \ref{fig.lpvgrid} illustrates the concept. A nonlinear system is linearized along a grid of 
Mach and altitude values, resulting in an array of linear systems. 
In this case the \texttt{pss} object is formed by combining a \texttt{ss} array of 
linearizations with an \texttt{rgrid} object representing the grid of parameter values 
(e.g. the grid of Mach and altitude values in Figure \ref{fig.lpvgrid}).
An alternative method of constructing grid-based LPV models in LPVTools 
is to use the \texttt{pgrid} atom (denoting a real scalar parameter defined on a grid of points).

<<LPVArray.png>>

\cite{Becker1993} or \cite{Wu1995}



The LPV system in Equation~\ref{eq:lpvsys} is
conceptually represented by a state-space system $S(\rho)$ that depends on a
time-varying parameter vector $\rho$ in some compact domain $\mathcal{P} \subset \R^{n_\rho}$.
A grid-based LPV model of this system is a collection of linearizations on a 
gridded domain of parameter values.
This is pictorially represented in Figure~\ref{fig.lpvgrid}, 
for an example system that depends on two parameters (in this case Mach number and altitude). 
For general LPV systems this conceptual representation requires storing the state
space system at an infinite number of points in the domain of $\rho$.
The grid based LPV approach, as implemented in LPVTools, approximates this conceptual
representation by storing the LPV system as a state space array
defined on a finite, gridded domain.
For each grid point $\hat{\rho}_k$ there is a corresponding LTI system 
$(A(\hat{\rho}_k),B(\hat{\rho}_k),C(\hat{\rho}_k),D(\hat{\rho}_k))$ which describes the dynamics of $S(\hat{\rho}_k)$ 
when $\hat{\rho}_k$ is held constant (note that $\hat{\rho}_k$, represents a constant vector corresponding to 
the $k$-th grid point, while $\rho_i$ is later used to denote the $i$-th element of the vector $\rho$.)
All the linearized systems on the grid have identical inputs $u$, outputs $y$ and state vectors $x$.
Together they form a LPV system approximation of $S(\rho)$.
This approach is motivated by the traditional gain-scheduling framework 
in aircraft flight control, for which models are typically constructed as 
linearizations around various flight operating conditions.

