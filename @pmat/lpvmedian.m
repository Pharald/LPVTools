function V = lpvmedian(M,PName)
% LPVMEDIAN  Median value along an independent variable of a PMAT object.
%
% V=LPVMEDIAN(M) computes the median value of M across all points in the 
% domain of M. V is a PMAT of the same dimension as M with V(i,j) 
% equal to the median of M(i,j) over all the independent variables in the 
% domain of M.
%
% V=LPVMEDIAN(M,PNAME) computes the median value of M along a particular
% independent variable specified in PNAME. PNAME must be a single 
% char array (string). V is a PMAT of the same dimension as M with V(i,j) 
% equal to the median of M(i,j) over the independent variable of M that was 
% specified in PNAME. 
%
% See also: median, lpvmax, lpvmin, lpvmean.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 2)

if nin==1
    % Median over all IVs
    if nout ==2
        error('No index returned for LPVMEDIAN over all IVs.');
    end

    % Reorder as [row col AD IV]
    niv = M.Domain.NumIV;
    nad = numel(size(M))-2;
    Mdata = permute(M.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);
    
    % Vectorize IV dimensions
    len = 2+nad+1;
    idx = cellstr(repmat(':',len,1))';       
    Mdata = Mdata(idx{:});
    
    % Median over IVs
    V = median(Mdata,len);
    V = pmat(V);    
else
    % Median over the single IV specified in PName
    
    % Check inputs
    if ~ischar(PName) && size(PName,1)==1
        error('Second input argument of LPVMEDIAN must be a single character string.');
    end
    
    % Find location of PName
    Dom = M.Domain;
    niv = Dom.NumIV;
    nad = numel(size(M))-2;
    idx=find(strcmp(PName,Dom.IVName));
    
    % Median over PName dimension
    Data = M.Data;
    dim = idx+2;
    VData=median(Data,dim);
    
    VData = permute(VData,[1:(dim-1) (dim+1):(2+niv+nad) dim]);
    
    % Remove PName from Domain
    Dom = lpvelimiv(Dom,PName);
    
    % Pack up result
    V = pmat(VData,Dom);
end

