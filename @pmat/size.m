function [varargout] = size(m,arg2)
% SIZE   Size of a PMAT object.
%
% For an M-by-N PMAT A with P independent variables and D array dimensions, 
% S = SIZE(A) returns a 1-by-(2+D) row vector. S(1) and S(2) are the row and
% column dimensions of A. For i>2, S(i) is the i-th array dimension of A.
% Use SIZE(A.Domain) to determine the dimensions of the P independent 
% variables.
%
% See also: size, length.

out = privatesize(m);
dsz = size(m.DomainPrivate);
[uidx,arrayidx] = upidx(m.DomainPrivate);
out = [out(1:2) dsz(arrayidx)];

if length(out)==3
   out = [out 1];
end

if nargin==1
    arg2 = nan;
end

nout = nargout;
if nout==2
    % Direct call to csize with nout = 2 put sthe product of the column
    % and IV dims into the second output. Instead SIZE with two 
    % output args should return only the row/col dimensions.
    tmpout = csize(out,arg2,nargin,3);
    varargout = {tmpout{1}, tmpout{2}};
else
    varargout = csize(out,arg2,nargin,nout);
end
