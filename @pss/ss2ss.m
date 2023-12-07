function sys = ss2ss(sys,T)
% SS2SS  Change of state coordinates for PSS objects
%
% SYS = SS2SS(SYS,T) performs the similarity transformation z=Tx on the
% state vector x of the PSS model.  The transformation T must be a 
% constant matrix.
%  
% see also: ss2ss, canon, balreal.

% NOTE PJS 3/24/2011: Revisit--Possibly implement PMAT T.

% Check # of input/output arguments
narginchk(2, 2)
if isa(T,'pmat')
    T = double(T);
    if ndims(T)>2
        error('Transformation T must be a constant matrix');
    end
end

try
    sys.DataPrivate = ss2ss(sys.DataPrivate,T);
catch ME
    throw(ME);
end
    
