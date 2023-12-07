%% LPVBALREAL - Gramian-based balancing for |pss| objects
%
%  
%% Syntax
%
%    SYSB = lpvbalreal(SYS) 
%    [SYSB,G] = lpvbalreal(SYS)
%    [SYSB,G,T,Ti] = lpvbalreal(SYS)
%    [SYSB,G,T,Ti] = lpvbalreal(SYS,...,INVERSE)
%
%% Description 
%
% |SYSB = lpvbalreal(SYS)| computes a balanced realization of the 
% parameter varying system |SYS|. |SYSB| is dervied by computing a single 
% balancing transformation for |SYS| and applying it at every point in its 
% domain. 
% 
% |[SYSB,G] = lpvbalreal(SYS)| returns |G|, a vector of singular values 
% describing the input-output mapping of |SYSB| (comparable to Hankel 
% singular values for LTI systems).
%
% |[SYSB,G,T,Ti] = lpvbalreal(SYS)| returnes the balancing transformation |T|, and
% its inverse |Ti|. 
%
% |[SYSB,G,T,Ti] = lpvbalreal(SYS,...,INVERSE)| provides an option to use a 
% alternative implementation of the algorithm which computes the balancing 
% transformation. If |INVERT| is True the alternative formulation is used. It
% can improve the accuracy of the solution for certain systems. The default 
% implementation assumes |INVERT=FALSE|. 