function V = lpvmean(M,PName)
% LPVMEAN  Mean value along an independent variable of a PMAT object.
%
% V=LPVMEAN(M) computes the mean value of M across all points in the domain 
% of M. V is a PMAT of the same dimension as M with  V(i,j) equal to the 
% mean of M(i,j) over all the independent variables in the domain of M.
%
% V=LPVMEAN(M,PNAME) computes the mean value of M along a particular
% independent variable specified in PNAME. PNAME must be a single 
% char array (string). V is a PMAT of the same dimension as M with V(i,j) 
% equal to the mean of M(i,j) over the independent variable of M that was 
% specified in PNAME. 
%
% See also: mean, lpvmax, lpvmin, lpvmean, lpvmedian

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 2)

if nin==1
    % Mean over all IVs
    if nout ==2
        error('No index returned for LPVMEAN over all IVs.');
    end

    % Reorder as [row col AD IV]
    niv = M.Domain.NumIV;
    nad = numel(size(M))-2;
    Mdata = permute(M.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);
    
    % Vectorize IV dimensions
    len = 2+nad+1;
    idx = cellstr(repmat(':',len,1))';       
    Mdata = Mdata(idx{:});
    
    % Mean over IVs
    V = mean(Mdata,len);
    V = pmat(V);    
else
    % Mean over the single IV specified in PName
    
    % Check inputs
    if ~ischar(PName) && size(PName,1)==1
        error('Second input argument of LPVMEAN must be a single character string.');
    end
    
    % Find location of PName
    Dom = M.Domain;
    niv = Dom.NumIV;
    nad = numel(size(M))-2;
    idx=find(strcmp(PName,Dom.IVName));
    
    % Mean over PName dimension
    Data = M.Data;
    dim = idx+2;
    VData=mean(Data,dim);
    
    VData = permute(VData,[1:(dim-1) (dim+1):(2+niv+nad) dim]);
    
    % Remove PName from Domain
    Dom = lpvelimiv(Dom,PName);
    
    % Pack up result
    V = pmat(VData,Dom);
end

