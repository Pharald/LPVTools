function a = depromote(a)
% DEPROMOTE   Demote the class of a PSTRUCT if possible.
%
% If a PSTRUCT S has no independent variables then B=depromote(S) will
% be a (standard) structure.
%
% See also: struct.

if isempty(a.Domain)
    a = a.DataPrivate;
end
