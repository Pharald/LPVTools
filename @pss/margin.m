function [Gm,Pm,Wcg,Wcp] = margin(SYS)
% MARGIN  Pointwise gain and phase margins, and crossover frequencies for a PSS
%  
% [Gm,Pm,Wcg,Wcp] = MARGIN(SYS) computes the gain margin Gm, the phase 
% margin Pm, and the associated frequencies Wcg and Wcp, for the SISO
% open-loop PSS SYS (continuous or discrete).  Gm, Pm, Wcg, Wcp are 
% returned as PMATs that provide the margins and crossover frequencies at
% each point in the domain of SYS. See DynamicSystem/MARGIN for the 
% definitions of the gain and phase margin.
%
% See also: margin.

% TODO PJS 5/28/2011: How do we handle allmargin and loopmargin?

Dom = SYS.Domain;
[Gm,Pm,Wcg,Wcp] = margin(SYS.Data);
Gm = pmat(shiftdim(Gm,-2),Dom);
Pm = pmat(shiftdim(Pm,-2),Dom);
Wcg = pmat(shiftdim(Wcg,-2),Dom);
Wcp = pmat(shiftdim(Wcp,-2),Dom);
