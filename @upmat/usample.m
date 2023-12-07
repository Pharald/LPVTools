function [B,Samples] = usample(UVars,N,varargin)
%USAMPLE   Generates random samples of UPMATs.
% 
%   B = USAMPLE(A,N) generates N random samples of the uncertainty in A 
%   and evaluates the UPMAT A at each one to generate a PMAT B. The output 
%   B is a [size(A) N] array of PMATs.
% 
%   B = USAMPLE(A) takes a single sample of the uncertainty and is
%   identical to B = USAMPLE(A,1).
%
%   [B,SAMPLES] = USAMPLE(A,N) returns the N samples of the uncerainty in A.
%   SAMPLES is a N-by-1 struct array, whose fields correspond to the
%   uncertainties in A. 
%
%   [B,SAMPLES] = USAMPLE(A,NAMES,N) only takes samples of the
%   uncertainties specified in NAMES. NAMES can be a string for a single 
%   uncertainty, or a cell array of strings for multiple uncertainties. If 
%   NAMES does not contain the names of all the uncertainties in A, the 
%   output B is a UPMAT containining the remaining uncertainties.
%
%   [B,SAMPLES] = USAMPLE(A,NAME1,N1,NAME2,N2,...) takes N1 samples of the
%   uncertainty NAME1, N2 samples of uncertainty NAME2, etc. The output B 
%   has size [size(A),N1,N2,...].
%
%   See also usample, usubs.

ni = nargin;
Dom = UVars.DomainPrivate;
Udata = UVars.DataPrivate;

if ni == 1
    N = 1;
end

if nargout==1
    B = usample(Udata,N,varargin{:});
else
    [B,Samples] =  usample(Udata,N,varargin{:});
end

if isuncertain(B)
	B = upmat(B,Dom);
else
    B = pmat(B,Dom);
end






    
