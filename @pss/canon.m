function [csys,T] = canon(sys,varargin)
% CANON  Pointwise canonical state-space realizations for PSS objects 
%
% CSYS = CANON(SYS,TYPE) computes a canonical state-space realization
% CSYS at each point in the domain of the PSS SYS. The string TYPE
% selects the type of realization and can be 'modal' (default) or
% 'companion'.
%
% [CSYS,T] = CANON(SYS,TYPE) also returns the transformation PMAT T that 
% relates the canonical state vector z to the original state vector
% x, by z = Tx at each point in the domain of SYS.
%
% CSYS = CANON(SYS,'modal',CONDT) specifies an upper bound CONDT on 
% the condition number of the block-diagonalizing transformation T. 
%
% See also: canon, ss2ss.

nout = nargout;
try
    csys = sys;
    npts = prod(sys.DomainPrivate.LIVData);        
    if nout==1
        for i=1:npts
            csys.DataPrivate(:,:,i) = canon(sys.DataPrivate(:,:,i),varargin{:});
        end
    else
        sz = [privatesize(sys) 1];
        ns = size(sys.DataPrivate.a,1);
        Data = zeros([ns ns sz(3:end)]);
        for i=1:npts
            [csys.DataPrivate(:,:,i),Data(:,:,i)] = ...
                canon(sys.DataPrivate(:,:,i),varargin{:});
        end
        T = pmat(Data,sys.DomainPrivate);
    end
catch ME
    throw(ME);
end
    
