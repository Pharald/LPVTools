function [Q,R,E] = qr(A,B)
% QR    Orthogonal-triangular decomposition for PMAT objects.
%
% [Q,R]=QR(A) computes the QR decomposition of A at each point in the 
% domain of A. See additional calling options in the QR help.
%
% See also: qr.

% AH 2/7/14 - See also: lu, null, orth, qrdelete, qrinsert, qrupdate?
% We don't have them implemented.

% NOTE PJS 3/23/2011: QR syntax for sparse matrices is not implemented.
% TODO PJS 4/2/2011: Implement functions listed in the "See also".

% Input error checking
nin = nargin;
nout = nargout;
narginchk(1, 2)
nargoutchk(0, 3)

szA = privatesize(A);
m = szA(1);
n = szA(2);
Data = A.DataPrivate;
if nin==1
    if nout <= 1
        X = zeros(szA); 
        for i=1:prod(szA(3:end))
            X(:,:,i) = qr(Data(:,:,i));
        end        
        Q = pmat(X,A.DomainPrivate);        
    elseif nout==2
        R = zeros([m n szA(3:end)]);
        Q = zeros([m m szA(3:end)]);
        for i=1:prod(szA(3:end))
            [Q(:,:,i),R(:,:,i)] = qr(Data(:,:,i));
        end        
        Q = pmat(Q,A.DomainPrivate);        
        R = pmat(R,A.DomainPrivate);         
    elseif nout==3
        R = zeros([m n szA(3:end)]);
        Q = zeros([m m szA(3:end)]);
        E = zeros([n n szA(3:end)]);
        for i=1:prod(szA(3:end))
            [Q(:,:,i),R(:,:,i),E(:,:,i)] = qr(Data(:,:,i));
        end        
        Q = pmat(Q,A.DomainPrivate);        
        R = pmat(R,A.DomainPrivate);         
        E = pmat(E,A.DomainPrivate);                 
    end 
        
elseif nin==2
    if isa(B,'pmat')
        B = double(B);
    end
    if ~isequal(B,0)
        error('Use qr(X,0) for economy size decomposition.')
    end
    if nout <= 1
        X = zeros(szA); 
        for i=1:prod(szA(3:end))
            X(:,:,i) = qr(Data(:,:,i),0);
        end        
        Q = pmat(X,A.DomainPrivate);        
    elseif nout==2
        if szA(1)<=szA(2)
            R = zeros([m n szA(3:end)]);
            Q = zeros([m m szA(3:end)]);
        else
            R = zeros([n n szA(3:end)]);
            Q = zeros([m n szA(3:end)]);
        end
        for i=1:prod(szA(3:end))
            [Q(:,:,i),R(:,:,i)] = qr(Data(:,:,i),0);
        end        
        Q = pmat(Q,A.DomainPrivate);        
        R = pmat(R,A.DomainPrivate);         
    elseif nout==3
        R = zeros([m n szA(3:end)]);
        Q = zeros([m m szA(3:end)]);
        E = zeros([n 1 szA(3:end)]);
        for i=1:prod(szA(3:end))
            [Q(:,:,i),R(:,:,i),E(:,:,i)] = qr(Data(:,:,i),0);
        end        
        Q = pmat(Q,A.DomainPrivate);        
        R = pmat(R,A.DomainPrivate);         
        E = pmat(E,A.DomainPrivate);                 
    end 
end

