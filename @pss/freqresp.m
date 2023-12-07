function [H,W] = freqresp(sys,varargin)
%  FREQRESP  Pointwise frequency response of a PSS.
%
% H = freqresp(SYS,...) computes the frequency response of the PSS SYS 
% at each point in the domain of SYS. H is a PMAT containing the frequency
% response of SYS at each point in the domain.
%
% See also: bode. 

narginchk(1,2)

if nargout == 1
    H = freqresp(sys.Data,varargin{:});
else
    error('Two output call not implemented')
end
niv = sys.Domain.NumIV;
nad = numel(size(sys.Data))-2-niv;
H = permute(H,[1 2 4:(3+niv) 3 (3+niv+1):(3+niv+nad)]);
H = pmat(H,sys.Domain);


