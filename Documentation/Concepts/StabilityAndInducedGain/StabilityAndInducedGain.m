%% Stability and Gain of an LPV system
% 
% LPVTools provides a suite of functions to analyze the stability and gain 
% of LPV systems. Meanwhile, LPVTools synthesis functions generate  
% controllers that are provide closed-loop stability for an LPV system, 
% while optimizing the gain.
% This section will discuss what stability and 
% gain mean for an LPV system. Furthermore, this section 
% highlights some of the computational issues that arise when LPV analysis
% conditions are implemented.
% 
%% Stability and Gain of an LPV system
%
% LPV systems are time-varying, state-space models of the form:
% 
% $$\left[ \begin{array}{c} \dot x (t) \\ y (t)\end{array} \right]
% = \left[ \begin{array}{cc} A(\rho(t)) & B(\rho(t)) \\ C(\rho(t)) & D(\rho(t))
% \end{array} \right] \left[ \begin{array}{c} x (t) \\ u (t)\end{array} \right]
% \ \ \ \ \ \ \ (1)$$
% 
% where $\rho \in \mathcal{R}^{n_\rho}$ is a vector of measurable parameters, 
% $y \in \mathcal{R}^{n_y}$ is a vector of outputs,
% $x \in \mathcal{R}^{n_x}$ is the state vector, $u \in \mathcal{R}^{n_u}$ is a vector 
% of inputs, and 
% $A\in \mathcal{R}^{n_x \times n_x}$, $B\in \mathcal{R}^{n_x \times n_u}$, 
% $C\in \mathcal{R}^{n_y \times n_x}$ and $D\in \mathcal{R}^{n_y \times n_u}$ are parameter 
% dependent matrices. 
% 
% The LPV system in Equation (1) depends on a set of time-varying parameters $\rho$. 
% The trajectories of the parameters are assumed to take on values in a 
% known compact set $\mathcal{P} \subseteq \mathcal{R}^{n_\rho}$, and to have known 
% bounds on their derivatives with respect to time: $\overline{\nu} \leq \dot{\rho} \leq \underline{\nu}$, 
% where $\overline{\nu}$ and $\underline{\nu} \in \mathcal{R}^{n_\rho}$.  
% A trajectory is said to be "rate unbounded" if $\overline{\nu} = \infty$ 
% and $\underline{\nu} = -\infty$.
% 
% The LPV system processes the inputs $u$ linearly, 
% but can depend nonlinearly on the time-varying parameter $\rho$. 
% The analysis problem is is to determine if the system is stable, and to 
% quantify the input-to-output gain of the system.
% Denote the LPV system in Equation (1) by $G(\rho)$.
% Analysis in the LPV framework determines if $G(\rho)$
% is internally exponentially stable, and whether the input/output map $G(\rho)$ from
% $u(t)$ to $y(t)$ has certain properties. 
%
%
% *Definitions of Gain*
% 
% LPVTools implements two methodologies for synthesis and analysis in the LPV framework. The two 
% methodologies differ in their formulation of the input/output map $G(\rho)$.
% The first methodology formulates this input/output map in
% terms of the induced $L_2$ norm (gain) of the system:
%
% $$ \| G(\rho) \|_{2 \to 2} = \max_{\rho \in \mathcal{P},~\overline{\nu}
% \leq \dot{\rho} \leq \underline{\nu}} 
% ~~\max_{u \in L_2,~\|u\|_2 \neq 0}
% \frac{\| G(\rho) u \|_2}{\| u \|_2} \ \ \ \ \ \ \ (2)
% $$
% 
% In calculating this induced norm it is assumed that $x(0)=0$.
% The second methodology formulates the input/output map in terms of the stochastic LPV bound 
% on $G(\rho)$:
% 
% $$
% stoch\left(G(\rho)\right) = \lim_{T\to\infty}  ~\max_{\rho \in \mathcal{P},~ 
% \overline{\nu} \leq \dot{\rho} \leq \underline{\nu}} ~
% E\left \lbrace \frac{1}{T}\int_{0}^{T} y^T(t)y(t) dt \right \rbrace  \ \ \ \ \ \ \  (3)
% $$
%
% which describes the variance of $y$ when the input $u$ is a zero mean, 
% white-noise processes with unit intensity. 


%% Computing the nominal $L_2$ norm of a grid-based LPV system:
% 
% |lpvnorm| implements algorithms to compute the gain of LPV
% systems. This section will review the analysis conditions that 
% |lpvnorm|  implements to compute the induced $L_2$ norm of a 
% grid-based nominal (not uncertain) LPV system. These analysis conditions
% will serve to illuminate many of the key issues in LPV analysis techniques.
% Refer to the references at the end of this chapter for conditions used in
% other analysis scenarios.
%
%
%%
% *The Objective*
% 
% The theory underpinning the LPV analysis results which are implemented 
% in |lpvnorm| frames the analysis problem in terms 
% of a  dissipation inequality. For the LPV system in Equation (1), 
% the problem boils down to a set Linear Matrix Inequalities (LMIs) 
% which need to be solved to prove that:
% 
% $$ \int_0^T y(t)^T y(t) \, dt \leq \gamma^2 \int_0^T u(t)^T u(t) \, dt \ \ \ \ \ \ \ (4)$$
% 
% for all $\rho \in \mathcal{P}$ and
% $\overline{\nu} \leq \dot{\rho} \leq \underline{\nu}$, with some 
% $\gamma \in \mathcal{R}^+$ and initial condition $x(0) = 0$.
% 
% Solving the LMIs to show that the dissipation inequality in Equation (4)
% holds, is sufficient to prove that the system is internally exponentially 
% stable, and that the gain of the system has a finite upper bound ($\gamma$). 
% The nominal induced $L_2$ norm analysis conditions used by |lpvnorm|
% are based on result by F. Wu. [1,2]

%%
% *Analysis Conditions*
% 
% The following theorem, taken from [1,2],
% gives a condition for an upper bound on the induced $L_2$ norm of the
% nominal LPV system $G(\rho)$ in Equation (1). For simplicity we will
% assume that the rate bounds on the parameter are symmetric: 
% $\nu = \overline{\nu} = -\underline{\nu}$.
% 
%%
% *_Theorem 1_*: If there exists a piecewise continuous symmetric function 
% $X:\mathcal{R}^{n_\rho} \rightarrow {\mathcal{R}^{n_x \times n_x}}$ and 
% a $\gamma \in \mathcal{R}^+$, such that $X(\rho)>0$  and
% 
% $$
% \left[ \begin{array}{ccc} A^T(\rho) X(\rho) + X(\rho) A(\rho)
% + \sum_{i=1}^{n_\rho} \beta_i \frac{\partial X}{\partial \rho_i} &
% X(\rho) B(\rho) & \gamma^{-1} C^T(\rho) \\
% B^T(\rho)X(\rho) & -I_{n_u} & \gamma^{-1} D^T(\rho) \\
% \gamma^{-1} C(\rho) &  \gamma^{-1} D(\rho) & -I_{n_y} \end{array} \right]<0
% \ \ \ \ \ \ \ (5)$$
% 
% $\forall \rho \in {\mathcal P}$, and $-\nu \leq \dot{\rho} \leq \nu$, 
% with $|\beta_i| \le \nu_i$ $(i=1,\ldots,n_\rho)$, then:
% 
% * The system $G$ is parametrically-dependent stable over ${\mathcal P}$.
% * $\exists k$ with $0\le k < \gamma$ such that $\|G\|_{2\to2} \le k$.
% 
% 
%%
% 
% The theorem above assume that the rate bounds of the time-varying parameter are symmetric, 
% but it can be extended to the unsymmetric case, and the software handles the unsymmetric case.
% The conditions in Theorem 1 are a parameterized set of
% linear matrix inequalities (LMIs) that must be verified for all $\rho
% \in {\mathcal P}$ and all $|\beta_i| \le \nu_i$. 
% The conditions are infinite dimensional, since $A(\rho)$, $B(\rho)$, 
% $C(\rho)$, $D(\rho)$ and $X(\rho)$ are all continuous functions of the 
% parameter $\rho$. 
% 
% *Implementation in LPVTools*
% 
% Its possible to obtain an approximate solution to the infinite 
% dimensional feasibility conditions in Theorem 1 by converting them into 
% a finite-dimensional set of Linear Matrix Inequalities (LMIs).
% This is accomplished by the following proceedure:
% 
% # Grid the set $\mathcal{P}$ into a set of $n_r$ 
% parameter values: $\{ \hat{\rho}_1, \hat{\rho}_2,...\hat{\rho}_{n_r}\}$.
% Require that the LMIs in Equation (5) hold at each grid point.
% # Pick a basis for $X(\rho)$ so that $X(\rho) = \sum_{k=1}^{n_b}f_k(\rho)X_k$, 
% where $n_b$ is the number of basis functions used to construct $X(\rho)$, 
% the scalar functions $f_1,\ldots, f_{n_b} : \mathcal{R}^{n_\rho} \to \mathcal{R}$ 
% are the chosen basis functions, and $X_1,\ldots,X_{n_b} \in \mathcal{R}^{n_x \times n_x}$ 
% are constant matrices to be determined 
% (see the <..\..\..\HowTo\BASISexample\html\BASISexample.html tutorial on picking basis functions> 
% for an example of how $f_1,\ldots, f_{n_b}$ are defined in LPVTools).
% If the parameter's in the LPV system are rate unbounded (i.e. $\nu =
% \infty$) then use a constant (parameter independent) Lyapunov matrix
% $X(\rho) = X \in \mathcal{R}^{n_x \times n_x}$.
% # Exploit the fact that the $\beta_i$
% enter affinely in Equation (4) to reduce the problem to $2^{n_\rho}$ LMIs at each grid point. 
% Specifically, if the LMIs hold for all combinations of 
% $\beta_i = \pm \nu_i$ (a total of $2^{n_\rho}$ combinations formed by 
% the $n_\rho$-dimensional polytope: $[-\nu_1,\nu_1] \times [-\nu_2,\nu_2] 
% \times \ldots \times [-\nu_{n_\rho},\nu_{n_\rho}]$) 
% then they hold for all $|\beta_i| \le \nu_i$. 
% This reduces the problem to $n_r 2^{n_\rho}$ LMIs total 
% ($n_r$ grid points, with $2^{n_\rho}$ LMIs at each point.) 
% # Solve for $\gamma$ and $X_1,\ldots,X_{n_b}$, subject to the $(n_r2^{n_\rho})$  
% LMIs formed at the grid points by the condition in Equation (5).
% 
% The function |lpvnorm| implements this proceedure to approximately solve 
% the conditions in Theorem 1 by enforcing the LMIs on the set of gridded 
% points in the domain of the grid-based LPV system (for a grid-based LPV 
% system the set of possible $\rho$ values, $\mathcal{P}$,
% is gridded as a matter of course during the modeling process). 
% 
%  
% The computational growth of these conditions is an issue.  Let $n_r$
% denote the total number of grid points used to approximate ${\mathcal P}$.
% A rate bounded analysis must enforce the LMI conditions at all $n_r$
% grid points and for all $2^{n_\rho}$ combinations of $\beta_i = \pm
% \nu_i$.  Thus there are a total of $n_r2^{n_\rho}$ constraints, each of
% dimension $(n_x+n_u+n_y)$.  If there are $n_b$ basis functions, then
% the Lyapunov matrix has $n_b$ symmetric matrix decision variables
% $\{X_j\}_{j=1}^{n_b}$ each of dimension $n_x \times n_x$. This gives a
% total of $n_r \frac{n_x(n_x+1)}{2}$ individual decision variables in
% the rate bounded analysis.  LMI optimization solvers have an
% asymptotic complexity that depends on both the number of decision
% variables and the number/dimension of LMI constraints.  For example,
% LMILab has a floating point operation growth of O($n_{row}n_v^3$) where
% $n_{row}$ is the total row dimension of the LMI conditions and $n_v$ is
% the total number of decision variables [3]. This
% complexity assumes the default Cholesky factorization of the Hessian
% matrix is used to solve the least squares problem that arises in each
% iteration.  Thus the complexity of solving the LPV analysis condition
% is roughly 
% $O\left( n_r2^{n_\rho}(n_x+n_u+n_y) \left(n_b n_x^2 \right)^3 \right)$.  
% This growth limits the analysis to a modest
% number of parameters, grid points, and basis functions.
% 
% *Alternative Approaches* 
%
% The LPV analysis problem is formulated differently when the system is 
% represented in the LFT-based LPV framework. In this case, the rate-bounds 
% can still be taken into account in the analysis, but they do not require 
% the user to define basis functions. The resulting feasability conditions 
% are different from the ones listed in the grid-based LPV analysis above.
% However, the implementations of the two approaches have many features in common: 
% Solution involves convex constraints (LMIs), 
% and the complexity grows with $O(2^{n_\rho})$.
% Further information on the analysis conditions for the LFT-based LPV approach 
% can be found in P. Apkarian and P.Gahinet [4], A. Packard [5], 
% A. Helmersson [6], and C. Scherer [7].
%
% The analysis conditions that apply for the stochastic LPV bound can be found in 
% the work by F. Wu [1], and the results for worst-case LPV analysis can be found
% in C. Scherer [7,8,9] and H. Pfifer and P. Seiler [10].


%% References
% 
% # F. Wu, "Control of Linear Parameter Varying Systems," PhD thesis, University of California,
% Berkeley, 1993.
% # F. Wu, X. Yang, A. Packard, and G. Becker, "Induced L2 norm control for LPV systems with
% bounded parameter variation rates," _International Journal of Nonlinear and Robust Control_,
% vol. 6, pp. 983-998, 1996.
% # P. Gahinet, A. Nemirovski, A. Laub, and M. Chilali, "LMI control toolbox user's guide," tech.
% rep., The Mathworks, 1995.
% # P. Apkarian and P.Gahinet, "A convex characterization of gain-scheduled Hinfinity controllers,"
% _IEEE Transactions on Automatic Control_, vol. 40, no. 5, pp. 853-864, 1995.
% # A. Packard, "Gain scheduling via linear fractional transformations," _Systems and Control
% Letters_, vol. 22, no. 2, pp. 79-92, 1994.
% # A. Helmersson, "An IQC-based stability criterion for systems with slowly
% varying parameters," Technical Report LiTH-ISYR-1979, Linkoping
% University 1997.
% # C. Scherer and S. Wieland, "Linear matrix inequalities in control,"
% Lecture notes for a course of the dutch institute of systems and
% control, Delft University of Technology, 2004.
% # C. Scherer and I. Kose, "Robustness with dynamic IQCs: An
% exact state-space characterization of nominal stability with
% applications to robust estimation," _Automatica_, Vol. 44, No. 7, 
% pp. 1666-1675, 2008.
% # C. Scherer, "LPV control and full-block multipliers," _Automatica_,
% Vol. 37, No. 3, pp. 361-375, 2001.
% # H. Pfifer, and P. Seiler. "Robustness analysis of linear parameter 
% varying systems using integral quadratic constraints," _International 
% Journal of Robust and Nonlinear Control_, 2014, doi: 10.1002/rnc.3240.

