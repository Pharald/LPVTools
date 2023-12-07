function [feas,Sfeas,xfeas] = iqcstab4(A,omega)
% Solve IQC LMI using cutting plane method
% The IQC LMI is
%      [G(jw); Im]'*S(jw)*[G(jw);Im] <= -tol*I
% where A = lft(Delta,G)
%
% INPUTS
%   A is an USS
%      (XXX Allow A to be a UFRD? FRDs with udyns are not allowed....)
%   omega is an Nw-by-1 array of frequency points.
%      (omega is only needed if A is a USS)
%
% OUTPUTS
%  feas: feas=1 if the LMI problem is feasible. Otherwise feas=0.
%  XXX: Output both G, S, and block in an info structure?
%
% MUSYNID

% Get dimensions
[ne,nd] = size(A);

% Call as a performance problem with a trivial perf. objective
% XXX Sfeas includes multiplier for zero performance block
[perfparm,Sfeas,xfeas] = iqcperf4(A,omega,'gain');
feas = ~isempty(perfparm);
    
