function out = lpvinterp(G,varargin)
% LPVINTERP  Interpolate a UPMAT
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,...) interpolates a UPMAT A on
% the domain specified by the NAME/VALUES pairs. Each VALUE is a vector of 
% values at which to interpolate A along the domain direction of the
% corresponding NAME. If an independent variable of A is not listed
% in the inputs, then all values along this domain direction are retained
% in the output B.
%
% B = LPVINTERP(A,NAME,VALUES) is an alternative syntax. NAME is an N-by-1
% cell array of characters and VALUES is an N-by-1 cell array of values.
%
% B = LPVINTERP(A,NAME1,VALUES1,NAME2,VALUES2,....,METHOD) includes a final 
% input argument called METHOD, which specifies the interpolation method to 
% be used. METHOD can be: 'nearest', 'linear', 'spline', or 'cubic'. The
% default is 'linear'. 
%
% See also: lpvsubs, lpvsplit, lpvsample.

% General TODO, akp.  Internal delays are in LFT form, and should be easy
% to handle - we just need access to the non-delay LFT part.  I/O delays
% are easy too, since they are just hanging off the ends, and should remain
% hanging of the end.   Ultimately, we need to work with TMW/Pascal to get this
% done properly.   We can probably write the correct method (IOModel) that
% effectively pulls out everything (InternalDelays, Delta, etc) and gives
% us an SS object, and a structure, and then we can pass back to an inverse
% function a modified version of ths SS and the structure, and it builds
% the USS back correctly.

% Get state-space data
try
   [M,Delta] = lftdata(G.Data);
catch
   error('Cannot interpolate matrices whose uncertainty description changes');
end
pM = pmat(M,G.Domain); 
out = lft(Delta,lpvinterp(pM,varargin{:}));
