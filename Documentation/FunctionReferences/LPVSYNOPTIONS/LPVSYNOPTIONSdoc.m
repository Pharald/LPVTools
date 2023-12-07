%% LPVSYNOPTIONS - Create a options object for LPV synthesis and analysis
%
%  
%% Syntax
%
%    opt = lpvsynoptions
%    opt = lpvsynOptions(Name1,Value1,Name2,Value2,...)
%
%% Description
% 
% |opt = lpvsynOptions(Name1,Value1,Name2,Value2,...)| creates a options
% object for parameter-varying synthesis and analysis. The |lpvsynOptions|
% object is used to specify the parameters of the optimization routines 
% used in the synthesis and analysis functions: |lpvsyn|, |lpvmixsyn|, 
% |lpvncfsyn|, |lpvloopshapesyn|, |lpvestsyn|, |lpvsfsyn|, |lpvnorm|,
% |lpvstochsyn|, |lpnvnorm|, |lpvwcgain|.
% 
% |opt = lpvsynOptions| creates an |lpvsynOptions| object initialized with 
% default values.
%
%
% The options are set using |NAME|, |VALUE| pairs, i.e. the options property 
% specified by the character string |NAME| is set to |VALUE|. 
% The setable options properties are specified below.  
% The default choice is specified in brackets.
% 
%--------------------------------------------------------------------------
%     NAME          |   VALUE     |    Description
%--------------------------------------------------------------------------
%   'Solver'        | ['lmilab']  | Optimization solver to be used.
%   -----------------------------------------------------------------------
%   'SolverOptions' | []          | Options passed directly to the solver.
%   -----------------------------------------------------------------------
%   'SolverInit'    | []          | Initial decision variables for LMI solver.
%   -----------------------------------------------------------------------
%   'Gammalb'       | [1e-6]      | Lower bound on closed-loop induced L2 norm.
%   -----------------------------------------------------------------------
%   'Gammaub'       | [1e6]       | Upper bound on closed-loop induced L2 norm.
%   -----------------------------------------------------------------------
%   'Xlb            | [1e-6]      | X Riccati variable lower bounded by Xlb*I
%   -----------------------------------------------------------------------
%   'Xub            | [1e6]       | X Riccati variable upper bounded by Xub*I
%   -----------------------------------------------------------------------
%   'Ylb            | [1e-6]      | Y Riccati variable lower bounded by Ylb*I
%   -----------------------------------------------------------------------
%   'Yub            | [1e6]       | Y Riccati variable upper bounded by Yub*I
%   -----------------------------------------------------------------------
%   'Method         | ['BackOff'] | String specifying the solution method:
%                   |             |----------------------------------------
%                   | 'MinGamma'  | Minimize the closed-loop induced L2 norm.
%                   |             |----------------------------------------
%                   | 'MaxFeas'   | Maximize the feasibility of the X Riccati 
%                   |             | subject to contraints Gammalb and Gammaub 
%                   |             | on the closed-loop induced L2 norm.
%                   |             |----------------------------------------      
%                   | 'BackOff'   | First solve the 'MinGamma' problem for 
%                   |             | GammaOpt and then solve a second stage
%                   |             | 'MaxFeas' problem with 
%                   |             | Gammaub = BackOffFactor*GammaOpt.
%                   |             | This two-stage solution improves
%                   |             | the numerical conditioning of the 
%                   |             | controller reconstruction.
%                   |             |----------------------------------------  
%                   | 'PoleCon'   | Constrain the closed-loop poles.
%   -----------------------------------------------------------------------              
%   'BackOffFactor' | [1.2]       | Multiplicative factor ( >= 1 ) to back off
%                   |             | the minimum gamma when Method = 'BackOff'.
%   -----------------------------------------------------------------------            

