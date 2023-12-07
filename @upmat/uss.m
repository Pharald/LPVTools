function out = uss(a,b,c,d,varargin)
% Call UPSS constructor
%
% SYS = uss(A,B,C,D) creates a uncertain parameter-varying state-space model 
%     on the domain of A, B, C and D. SYS will be a UPSS.
%
% SYS = uss(A) creates a parameter dependent static gain matrix on the 
%     domain of A.
% 
% See also: uss, upss.

if nargin==1
    out = upss(a);
elseif nargin==4
    out = upss(a,b,c,d);
else
    out = upss(a,b,c,d,varargin{:});
end
