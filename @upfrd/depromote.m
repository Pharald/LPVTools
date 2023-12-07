function a = depromote(a)
% DEPROMOTE   Demote the class of a UPFRD if possible.
%
% If a UPFRD SYS has no independent variables then B=depromote(SYS) will
% be an uncertain frequency response.
%
% See also: ufrd.

if isempty(a.Domain)
    a = a.DataPrivate;
end
