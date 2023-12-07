function a = depromote(a)
% DEPROMOTE   Demote the class of a UPMAT if possible.
%
% If a UPMAT SYS has no independent variables then B=depromote(SYS) will
% be an uncertain frequency response.
%
% See also: umat.

if isempty(a.Domain)
    a = a.DataPrivate;
end
