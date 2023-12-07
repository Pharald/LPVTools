function out = order(m)
% ORDER  Compute the system order for UPSS
%
% N = ORDER(SYS) returns the order of the system SYS. For a
%      UPSS, N corresponds to the number of states in SYS.
%
% See also: order.

out = numel(m.Data.StateName);

end