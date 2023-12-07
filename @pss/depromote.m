function a = depromote(a)
% DEPROMOTE   Demote the class of a PSS if possible.
%
% If a PSS SYS has no independent variables then B=depromote(SYS) will
% be a state-space system.
%
% See also: ss.

if isempty(a.Domain);
    a = a.DataPrivate;
end
