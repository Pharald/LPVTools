function sysout = tmpresid(sys,poleloc,method)
% TMPRESID   Remove high frequency poles from a PLFTSS.
% 
% P = TMPRESID(G,ELIM) eliminates the number of high frequency poles
% associated with ELIM.

% function POUT = tmpresid(P,ELIM)
%
% POUT = tmpresid(P,ELIM) 
%
% POUT = tmpresid(SYS,ELIM,METHOD) also specifies the state elimination
% method. Available choices for METHOD include 'MatchDC' (default) or
% 'Truncate'.
%


% **TO BE REPLACED**
% This is a temporary model reduction algorithm to implement truncation and 
% residualization for PLFTSS objects. Need to truncate or residualized high 
% frequency modes of a PLFTSS object much like MODRED. 
%
% This is a bit tricky as the way to describe the parameter dependendence
% of a system in LFT form, is to write it out as the LFT of a parameter 
% block, and an integrator block,
%
%  x_dot = -rho*x +u, y = x, 
%
% which can be expressed as the upper-LFT 
%
% M = [0,1,1; 1,0,0; 1,0,0]
% Delta = [rho,0; 0 1/s]
%
% Hence, poles of the system are dependent on the parameter, and different 
% parameter values yield different results. 
%
% NEEDS TO BE REPLACED

% 6 Feb 14: GJB 

if nargin<=3
    method =  'MatchDC';
end

idx = logical([ones(poleloc,1); zeros(11,1)]); 
%parms = sys.Parameter;
%parmname = fieldnames(parms);
[MDelta,RHO] = lftdata(sys,[],'Parameters');

MDR = canon(MDelta,'modal'); % modal form, but max ev at the top
MDR = modred(MDR,idx,method);
sysout = lft(RHO,MDR);
