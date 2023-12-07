function [varargout] = size(m,arg2)
% SIZE   Size of a PSS object.
%
% For an M-by-N PSS SYS with P independent variables and D array dimensions, 
% S = SIZE(SYS) returns a 1-by-(2+D) row vector. S(1) and S(2) are the row 
% and column dimensions of SYS. For i>2, S(i) is the i-th array dimension 
% of SYS. Use SIZE(SYS.Domain) to determine the dimensions of the P 
% independent variables.
%
% See also: size, length.

%JT 15.07.2016 added compatability for size(m,'order') calls
if nargin==2 && isa(arg2,'char') && strcmp(arg2,'order')  
    varargout{1} = order(m);
    return
end

out = privatesize(m);
dsz = size(m.DomainPrivate);
[uidx,arrayidx] = upidx(m.DomainPrivate);
out = [out(1:2) dsz(arrayidx)];
% TODO: AKP there has to be an easier way to properly strip the trailing 1s
if length(out)==3
   if out(3)~=1
      out = [out 1];
   else
      out(3) = [];
   end
end
% strip off trailing 1's, with rule that size is 2 or >=4
while (numel(out)==3 || numel(out)>4) && out(end)==1
   out(end) = [];
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

