function out = order(m)
% ORDER  Compute the system order for PSS
%
% N = ORDER(SYS) returns the order of the system SYS. For a
%      PSS, N corresponds to the number of states in SYS.
%
% See also: order.

out = numel(m.Data.StateName);

end