%% Model Reduction for LPV systems
% A primer on model reduction in the LPV framework.

%% Introduction
% 
% LPVTools provides tools for LPV model reduction. LPV model reduction is
% different from Linear Time-Invariant (LTI) model reduction techniques 
% which act on a single point, because they perform the model reduction for
% all values of the scheduling parameter simultaneously. The resulting
% reduced order model is still a LPV model with consistent state, input and
% output vectors. If LTI model reduction techniques (e.g. |balreal|) are
% applied to a LPV model, the resulting model may lose 
% <..\..\Concepts\StateConsistency\html\StateConsistency.html state consistency>
% and the resulting reduced order model is no longer a LPV system.
% LPVTools provides two functions for LPV model reduction. |lpvbalancmr|
% performs balanced trunctation, and provides the option of weighting 
% different frequency bands in the model reduction to emphasize accuracy 
% for some dynamics while de-emphasizing others.. However it is 
% restricted to stable LPV systems. |lpvncfmr| performs a contractive 
% coprime factorization of a LPV system, and can handle unstable LPV
% systems.
% 
% 
% *Further Reading*
% 
% # G. D. Wood, "Control of parameter-dependent mechanical systems," Ph.D. 
% Dissertation, University of Cambridge, 1995.
% #  G. D. Wood, P. J. Goddard, and K. Glover, "Approximation of linear 
% parameter-varying systems," _IEEE Conference on Decision and Control_,
% Vol. 1, pp 406-411, 1996.
% # R. Widowati, R. Bambang, R. Sagari, S. M. and Nababan, 
% “Model reduction for unstable LPV system based on coprime
% factorizations and singular perturbation,” _5th Asian Control
% Conference_, Vol. 2, pp. 963-970, Melbourne, Australia, 2004.

%% LPV Model Reduction Commands
% 
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVGRAM\html\LPVGRAMdoc.html">LPVGRAM</a> </td>
% <td>Compute Gramians for PSS objects.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVBALREAL\html\LPVBALREALdoc.html">LPVBALREAL</a> </td>
% <td>Perform Gramian-based balancing for PSS objects.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVBALANCMR\html\LPVBALANCMRdoc.html">LPVBALANCMR</a> </td>
% <td>Balanced truncation model reduction.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVNCFMR\html\LPVNCFMRdoc.html">LPVNCFMR</a></td>
% <td>Balanced truncation model reduction through contractive
%   coprime factorization.</td>
% </tr>
% </table>
% </html>

%% Examples and How To
% 
% * <matlab:open(fullfile(docroot,'robust/simplify-models.html')) LTI Model Reduction>
% * <..\..\HowTo\StableModelReduction\html\StableModelReduction.html Model Reduction for a stable LPV system>
% * <..\..\HowTo\UnstableModelReduction\html\UnstableModelReduction.html Model Reduction for an unstable LPV system>

%% Concepts
% 
% * <..\..\Concepts\StateConsistency\html\StateConsistency.html State
% Consistency>
