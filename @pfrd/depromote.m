function a = depromote(a)
% DEPROMOTE   Demote the class of a PFRD if possible.
%
% If a PFRD SYS has no independent variables then B=depromote(SYS) will
% be a frequency response.
%
% See also: frd.

if isempty(a.Domain)
    a = a.DataPrivate;
end
