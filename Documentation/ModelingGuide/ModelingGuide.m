%% Modeling LPV systems
% A primer on modeling LPV systems.

%% LPV Modeling Commands
% 
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\ClassReference\PGRID\html\PGRIDdoc.html">PGRID</a> </td>
% <td>Gridded real parameter.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\RGRID\html\RGRIDdoc.html">RGRID</a> </td>
% <td>Rectangular grid of parameter values.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\PMAT\html\PAMTdoc.html">PMAT</a> </td>
% <td>Parameter-varying matrix.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\PSS\html\PSSdoc.html">PSS</a> </td>
% <td>Parameter-varying state-space system.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\PFRD\html\PFRDdoc.html">PFRD</a> </td>
% <td>Parameter-varying frequency response data model.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\UPMAT\html\UPMATdoc.html">UPMAT</a> </td>
% <td>Parameter-varying uncertain matrix.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\UPSS\html\UPSSdoc.html">UPSS</a> </td>
% <td>Parameter-varying uncertain state-space system.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\UPFRD\html\UPFRDdoc.html">UPFRD</a> </td>
% <td>Parameter-varying uncertain frequency response data model.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\BASIS\html\BASISdoc.html">BASIS</a> </td>
% <td>Parameter-varying basis function for analysis and synthesis.</td>
% </tr>
% <tr>
% <td><a href="..\..\ClassReference\PSTRUCT\html\PSTRUCTdoc.html">PSTRUCT</a> </td>
% <td>Parameter-varying structure.</td>
% </tr>
% </table>
% </html>
%
% <html>
% <table border=1>
% <tr>
% <td><a href="..\..\FunctionReferences\DOMUNION\html\DOMUNIONdoc.html">DOMUNION</a> </td>
% <td>Map LPV objects onto a common domain.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSPLIT\html\LPVSPLITdoc.html">LPVSPLIT</a> </td>
% <td>Extract LPV model data from a subset of its parameter domain.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVINTERP\html\LPVINTERPdoc.html">LPVINTERP</a> </td>
% <td>Interpolate a grid-based LPV model.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSUBS\html\LPVSUBSdoc.html">LPVSUBS</a> </td>
% <td>Substitute values of parameters.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVELIMIV\html\LPVELIMIVdoc.html">LPVELIMIV</a> </td>
% <td>Eliminate parameters which only have a single grid point.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVSAMPLE\html\LPVSAMPLEdoc.html">LPVSAMPLE</a> </td>
% <td>Sample a grid-based LPV object.</td>
% </tr>
% <tr>
% <td><a href="..\..\FunctionReferences\LPVBALANCE\html\LPVBALANCEdoc.html">LPVBALANCE</a> </td>
% <td>Diagonal scaling for LPV models.</td>
% </tr>
% </table>
% </html>


%% Examples and How To
% 
% * <..\..\GettingStarted\MiniTutorials\Grid_Modeling\html\Grid_Modeling.html Tutorial: Constructing grid-based LPV models>
% * <..\..\GettingStarted\MiniTutorials\LFT_Tutorial\html\LFT_Tutorial.html Tutorial: Constructing LFT-based LPV models>
% * <..\..\GettingStarted\MiniTutorials\GridToLFT\html\GridToLFT.html Tutorial: Conversion between LFT and LPV models>
% * <..\..\HowTo\LPV_fromAnalyticalLinearization\html\LPV_fromAnalyticalLinearization.html Tutorial: Creating a grid-based LPV model from analytical linearization.>
% * <..\..\HowTo\LPV_fromNonlinearSim\html\LPV_fromNonlinearSim.html Tutorial: Creating a grid-based LPV model from a nonlinear model.>
 

%% Concepts
% 
% * <..\..\Concepts\PermissibleTrajectories\html\PermissibleTrajectories.html Permissible
% Parameter Trajectories>
% * <..\..\GettingStarted\LPVSystems\html\LPVSystems.html#2 Grid-based LPV model.>
% * <..\..\GettingStarted\LPVSystems\html\LPVSystems.html#3 LFT-based LPV model.>
% * <..\..\Concepts\QuasiLPV\html\QuasiLPV.html Quasi-LPV Models>
% * <..\..\Concepts\StateConsistency\html\StateConsistency.html State
% Consistency>
