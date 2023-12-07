function a = depromote(a)
% DEPROMOTE   Demote the class of a PMAT if possible.
%
% If a PMAT M has no independent variables then B=depromote(M) will
% be a double matrix.
%
% See also: double.

if isempty(a.Domain)
    a = a.DataPrivate;
end
