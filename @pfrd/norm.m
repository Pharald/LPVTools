function [n,f] = norm(sys,P,tol)
% NORM   Dynamic system norm of a PFRD object.
%
% NORM(S) is the H2 norm of S at each point in the domain of S. NORM(S,2)
% also returns the H2 norm.
%
% NORM(S,inf) is the Hinfty norm of S at each point in the domain of S. 
% NORM(S,inf,TOL) specifies the relative tolerance for the computed norm.
%
% [NINF,FPEAK] = NORM(S,inf) returns PMATs NINF and FPEAK.  NINF if the
% Hinfty norm and FPEAK the frequency at which the system achives the peak 
% for each point in the domain of S.
%
% See also: norm, sigma.

% TODO PJS 4/29/2011: Allow P and/or TOL to be a PMAT? The SS/NORM allows
% SS arrays for inputs and thus this PSS/NORM code does not need a  
% for-loop. However, SS/NORM does not allow P or TOL to be an array. Thus 
% allowing P/TOL to be varying in this PSS/NORM would require a for-loop.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 3)
nargoutchk(0, 2)
if nin>=2 && isa(P,'pmat');
    P = double(P);
    if ndims(P)>2
        error('Varying P is not allowed.');
    end
end
if nin==3 && isa(tol,'pmat');
    tol = double(tol);
    if ndims(tol)>2
        error('Varying TOL is not allowed.');
    end
end

szm = privatesize(sys);
if nin==1
    n = norm(sys.DataPrivate);    
    n = pmat( reshape(n,[1 1 szm(3:end)]) , sys.DomainPrivate);
    f = pmat;
elseif nin==2
    if nout==1        
        n = norm(sys.DataPrivate,P);
    else
        [n,f] = norm(sys.DataPrivate,P);
        f = pmat( reshape(f,[1 1 szm(3:end)]) , sys.DomainPrivate);
    end
    n = pmat( reshape(n,[1 1 szm(3:end)]) , sys.DomainPrivate);
else    
    if nout==1        
        n = norm(sys.DataPrivate,P,tol);
    else
        [n,f] = norm(sys.DataPrivate,P);
        f = pmat( reshape(f,[1 1 szm(3:end)]) , sys.DomainPrivate);
    end
    n = pmat( reshape(n,[1 1 szm(3:end)]) , sys.DomainPrivate);
end

