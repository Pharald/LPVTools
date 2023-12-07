function rsys = modred(sys,elim,method)
% MODRED  Model reduction for PSS objects
%
% RSYS = MODRED(SYS,ELIM) reduces the order of the PSS object
% by eliminating the states specified by the vector ELIM. The reduction
% is done at each point in the domain of SYS.
%  
% RSYS = MODRED(SYS,ELIM,METHOD) specifies the method used for the
% state elimination. METHOD can be 'MatchDC' or 'Truncate'.
%
% See also: modred, balreal.

try
    rsys = sys;
    npts = prod(sys.DomainPrivate.LIVData);        
    if nargin==2
        for i=1:npts
            rsys.DataPrivate(:,:,i) = modred(sys.DataPrivate(:,:,i),elim);
        end
    else
        for i=1:npts
            rsys.DataPrivate(:,:,i) = modred(sys.DataPrivate(:,:,i),elim,method);
        end
    end
catch ME
    throw(ME);
end
    
