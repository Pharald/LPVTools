function C = binop(A,B,op,varargin)
% BINOP Binary operations for a PSTRUCT
% 
% BINOP is a utility function that handles binary operations (e.g. blkdiag
% horzcat, etc.) on PSTRUCTs. The user should not call BINOP directly. 
% Instead each binary operation has a front end function that the user 
% calls (e.g. horzcat). The frontend function then call calls BINOP to do 
% the required operations. 

%% PSTRUCT/BINOP

% NOTE PJS 4/7/2011: FEVAL is a bit slower than directly performing the 
% operation, eg.
%   if M is a matrix then M*M is about 4x faster than
%   feval('mtimes',M,M) and M*M is about 2x faster than
%   fh(M,M) where fh=@mtimes.
% The main BINOPs are implemented directly (e.g. plus calls double/+).


% Use switchyard to handle lifting if classes are different
if ~isequal(class(A),class(B))
   [A,B]=switchyard(A,B);
   nin = nargin;
   if nin==3
      C=binop(A,B,op);
   else
      C=binop(A,B,op,varargin{:});
   end
   return
end

% Define A and B on a common domain 
% (domunion will lift to pmat if required)
[Aext,Bext] = domunion(A,B);
szA = [privatesize(Aext) 1];
Adata = Aext.DataPrivate;
Bdata = Bext.DataPrivate;


% Perform binary operation at each point in the combined domain
switch op
    case 'horzcat'
        % TODO PJS 6/28/2011: [A,[]] will error if A is a nontrivial PMAT.
        Cdata = cat(2,Adata,Bdata);
    case 'vertcat'
        Cdata = cat(1,Adata,Bdata);
    
    case 'blkdiag'
        szB = [privatesize(Bext) 1];
        Z12 = zeros([szA(1),szB(2:end)]);
        Z21 = zeros([szB(1),szA(2:end)]);
        Cdata = [Adata Z12; Z21 Bdata];    
    otherwise
        % TODO - Fix this
        error('This binop is not iplimented yet for PSTRUCTs')
end
C = pstruct(Cdata,Aext.DomainPrivate);
