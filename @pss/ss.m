function out = ss(a,b,c,d,varargin)
% Call PSS constructor
%
% SYS = ss(A,B,C,D) creates a parameter-varying state-space model on  
%     the domain of A, B, C and D. SYS will be a PSS.
%
% SYS = ss(A) promotes A to a PSS if its a SS.
% 
% See also: ss, pss.

if nargin==1
    out = pss(a);
elseif nargin==4
    out = pss(a,b,c,d);
else
    out = pss(a,b,c,d,varargin{:});
end