# LPVTools

## Description
The Linear Parameter-Varying Toolbox (LPVTools) is a toolbox for modeling, analysis and synthesis in the Linear Parameter-Varying framework (LPV). LPV framework provides a mathematically rigorous approach to the design of gain-scheduled controllers, with powerful guarantees on their performance and robustness in the presence of time-varying dynamics and uncertainty.

LPVTools provides LPV data structures and a set of tools for modeling, simulation, analysis and synthesis in the LPV framework. Its capabilities include tools for synthesis of parameter-varying output feedback controllers, state-feedback controllers, and estimators, which yield optimized performance for a given set of allowable parameter trajectories. Tools are provided for analysis of the stability and input-to-output gain of LPV systems (with and without uncertainty). Tools are provided for performing model reduction on LPV models. And finally, tools are provided for simulating the time-domain response of LPV systems along user-supplied parameter trajectories.

## Key Features
- Modeling of parameter dependent systems and gain-scheduled control laws.
- LPV analysis and control synthesis.
- Simulation of LPV systems.
- Model reduction for parameter-dependent systems.

## System Requirements
LPVTools requires MATLAB®, Simulink®, the Control System Toolbox®, and the Robust Control Toolbox®. LPVTools makes use of the Control System and Robust Control Toolbox’s data structures, control synthesis and analysis algorithms.

## Installation
To add LPVTools to MATLAB path run the addlpv script

## Support
For any suggestions on future enhancement or nice to have features of these tools and reporting bugs, please email Harald Pfifer (harald.pfifer@tu-dresden.de) or Emily Burgin (emily.burgin@tu-dresden.de)

## Authors and acknowledgment

Current Developers: Harald Pfifer, Emily Burgin

The toolbox has been originally developed by G. Balas, A. Hjartarson, A. Packard, and P. Seiler. We are grateful for them to allow us taking over working on these tools.

## License
GNU Affero General Public License v3.0

## Caution
The user is cautioned that computer programs developed as a part of this research may not have been exercised for all cases of interest. While every effort is being made, within available time, to ensure that the programs are free of computational and logical errors, they cannot be considered validated. Any application of these programs without additional verification is at the risk of the user.
