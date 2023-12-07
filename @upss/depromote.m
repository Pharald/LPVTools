function a = depromote(a)
% DEPROMOTE   Demote the class of a UPSS if possible.
%
% If a UPSS SYS has no independent variables then B=depromote(SYS) will
% be an uncertain state-space system.
%
% See also: uss.

if isempty(a.Domain)
    a = a.DataPrivate;
end
