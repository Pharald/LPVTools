% LPVTools Toolbox
% Version 1.00 (R2014b) 21-April-2015
% 
%
% Grid based parameter dependent matrices and systems
%   pgrid         - Grid of parameter values
%   rgrid         - Rectangular grid of parameter values.
%   pmat          - Parameter-varying dependent matrix.
%   pss           - Parameter-varying state-space model.
%   pfrd          - Parameter-varying frequency response data model.
%   upmat         - Parameter-varying uncertain matrix.
%   upss          - Parameter-varying uncertain state-space model.
%   upfrd         - Parameter-varying uncertain frequency response data model.
%   basis         - Basis function for rate-bounded analysis and synthesis.
%   pstruct       - Parameter dependent structure.
%
% LFT based parameter dependent matrices and systems
%   tvreal        - Time-varying real parameter.
%   plftmat       - Parameter-varying matrix in LFT framework.
%   plftss        - Parameter-varying ss array in LFT framework.
%
% Manipulation of parameter dependent models
%   Domunion      - Map grid-based LPV objects onto a common domain
%   grid2lft      - Transform a grid-based LPV model into a LFT model.
%   lft2grid      - Transform a LFT model into a grid-based LPV model.
%   lpvsplit      - Extract grid-based LPV model data based on IV range.
%   lpvsubs       - Substitutes values for parameters.
%   lpvelimiv     - Eliminate singleton independent variables.
%   lpvsample     - Sample a grid-based LPV object.
%   lpvinterp     - Interpolate value of LPV object between grid points.
%   lpvbalance    - Diagonal scaling for PSS objects.
%   lpvplot       - Plot parameter dependent matrix data over IV range.
%
% Model order reduction
%   lpvgram       - Compute Gramians for PSS objects.
%   lpvbalreal    - Perform Gramian-based balancing for PSS objects.
%   lpvbalancmr   - Balanced truncation model reduction for PSS objects.
%   lpvncfmr      - Balanced normalized coprime factor model reduction.
%   
% Robustness and worst-case analysis.
%   lpvnorm       - Bound on induced L2 norm for PSS systems.
%   lpvwcgain     - Worst-case gain of an uncertain LPV system.
%
% Controller synthesis.
%   lpvsyn        - Synthesize an LPV controller.
%   lpvstochsyn   - Synthesize a LPV controller for a stochastic LPV system.
%   lpvncfsyn     - Normalized coprime factor LPV controller synthesis.
%   lpvsfsyn      - Synthesize a LPV state feedback controller.
%   lpvestsyn     - Synthesize a LPV state estimator.
%   lpvmixsyn     - Parameter-varying loop-shaping synthesis.
%   lpvloopshape  - Parameter-varying mixed-sensitivity synthesis.
%   lpvsynOptions - Create a options object for LPV synthesis and analysis
%
% Time-domain analysis
%   lpvlsim       - Simulate parameter dependent time response of a PSS.
%   lpvstep       - Simulate parameter dependent step response of a PSS.
%   lpvinitial    - Simulate initial conditions response of a PSS.
%   lpvimpulse    - Simulate parameter dependent impulse response of a PSS.
%
% Simulink
%   LPVBlockLibrary - Simulink block library for LPV models.
%
% Documentation.
%   Type "doc" and choose Supplemental Software for access to user manual.

