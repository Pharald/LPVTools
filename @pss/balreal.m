function [sysb,G,T,Ti] = balreal(sys,varargin)
% BALREAL   Pointwise Gramian-based balancing for PSS objects
%
% [SYSB,G] = BALREAL(SYS) computes a balanced state-space realization
% SYSB for the stable portion of the PSS SYS. This balancing is done at
% each point in the domain of SYS. G is a PMAT containing the vector of
% Hankel Singular Values of SYS at each point in the domain.
%
% [SYSB,G,T,Ti] = BALREAL(SYS,OPTIONS) allows options to be set.
% See HSVDOPTIONS for details on balancing options.  In addition,
% T and Ti are the balancing similarity transformation and its inverse at
% each point in the domain of SYS.
%
% See also: balreal, lpvbalreal, hsvdOptions, gram, modred.

nout = nargout;
try
    sysb = sys;
    Domain = sys.DomainPrivate;
    pts = Domain.LIVData;
    pts = pts(:)';
    npts = prod(pts);
    ns = size(sys.DataPrivate.a,1);
    if nout==1
        for i=1:npts
            sysb.DataPrivate(:,:,i) = balreal(sys.DataPrivate(:,:,i),varargin{:});
        end
    elseif nout==2
        G = zeros(ns,1,npts);
        for i=1:npts
            [sysb.DataPrivate(:,:,i), G(:,:,i)] = ...
                balreal(sys.DataPrivate(:,:,i),varargin{:});
        end
        G = pmat( reshape(G,[ns 1 pts]), Domain);
    else
        G = zeros(ns,1,npts);
        T = zeros(ns,ns,npts);
        Ti = zeros(ns,ns,npts);
        for i=1:npts
            [sysb.DataPrivate(:,:,i), G(:,:,i), T(:,:,i), Ti(:,:,i)] = ...
                balreal(sys.DataPrivate(:,:,i),varargin{:});
        end
        G = pmat( reshape(G,[ns 1 pts]), Domain);
        T = pmat( reshape(T,[ns ns pts]), Domain);
        Ti = pmat( reshape(Ti,[ns ns pts]), Domain);
    end
catch ME
    throw(ME);
end
