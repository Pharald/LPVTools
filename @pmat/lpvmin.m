function [V,I] = lpvmin(M,PName)
% LPVMIN  Smallest component along independent variables of a PMAT object.
%
% V=LPVMIN(M) computes the smallest element of M across all points in the 
% domain of M.  V is returned as a PMAT with no independent variables.
% V is a PMAT of the same dimension as M with V(i,j) equal to the minimum 
% of M(i,j) over all the independent variables in the domain of M.
%
% V=LPVMIN(M,PName) computes the smallest element of M along the particular
% independent variable specified in PName. PNAME must be a single char
% array (string). V is returned as a PMAT with the independent variable 
% specified in PNAME removed. V is a PMAT of the same dimension 
% as M with V(i,j) equal to the minimum of M(i,j) over the independent 
% variable of M that was specified in PNAME. 
%
% [V,I]=LPVMIN(M,PName) returns an additonal output I, which is a PMAT of  
% the same dimension as M containing the indices (single-index) of the  
% minimum values. The indices refer to the dimension specified in PName.
%
% See also: min, max, lpvmax, lpvmean, lpvmedian.

% Check # of input/output arguments
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 2)

if nin==1
    % Min over all IVs
    if nout ==2
        error('No index returned for LPVMIN over all IVs.');
    end

    % Reorder as [row col AD IV]
    niv = M.Domain.NumIV;
    nad = numel(size(M))-2;
    Mdata = permute(M.Data,[1 2 (3+niv:2+niv+nad) (3:2+niv)]);
    
    % Vectorize IV dimensions
    len = 2+nad+1;
    idx = cellstr(repmat(':',len,1))';       
    Mdata = Mdata(idx{:});
    
    % Min over IVs
    V = min(Mdata,[],len);
    V = pmat(V);    
else
    % Min over the single IV specified in PName
    
    % Check inputs
    if ~ischar(PName) && size(PName,1)==1
        error('Second input argument of LPVMIN must be a single character string.');
    end
    
    % Find location of PName
    Dom = M.Domain;
    niv = Dom.NumIV;
    nad = numel(size(M))-2;
    idx=find(strcmp(PName,Dom.IVName));
    
    % Min over PName dimension
    Data = M.Data;
    dim = idx+2;
    [VData,IData]=min(Data,[],dim);
    
    VData = permute(VData,[1:(dim-1) (dim+1):(2+niv+nad) dim]);
    IData = permute(IData,[1:(dim-1) (dim+1):(2+niv+nad) dim]);
    
    % Remove PName from Domain
    Dom = lpvelimiv(Dom,PName);
    
    % Pack up result
    V = pmat(VData,Dom);
    I = pmat(IData,Dom);
end

