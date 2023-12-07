function [dMd,dr,dci] = lpvbalance(sys,blk)
% LPVBALANCE   Diagonal scaling for PSS objects
%
% [P,D] = lpvbalance(S) computes a single diagonal similarity transformation 
% to improve the conditioning of the Ny-by-Nu PSS S at all points in the
% domain. The transformation D is computed using a generalized version of
% Osborne's iteration. D is returned as an Nx-by-Nx diagonal, double 
% matrix where Nx is the state dimension of S. The scaled PSS P is 
% obtained by applying the similarity transformation D at each point in 
% the domain of S.
%
% [P,DR,DC] = lpvbalance(S,BLK) computes diagonal transformations applied to
% the states and input/output channels. The state matrices of the scaled 
% PSS P are obtained from DR*[A B; C D]*DC where A, B, C, D are the state 
% matrices of S. BLK is a K-by-2 matrix that specifies the partitioning 
% dimensions of DR and DC. If BLK = [c1 r1; ... ; ck rk] then DR and DC 
% are partitioned as:
%    DR = blkdiag( D, d1*I_r1, ...,  dk*I_rk )
%    DC = blkdiag( inv(D), (1/d1)*I_c1, ...,  (1/dk)*I_ck )
% where D is a diagonal, Nx-by-Nx matrix, and the notation I_r1= eye(r1), 
% represents the r1-by-r1 identity matrix. The block partitioning must be 
% consistent with the input/output dimensions of S, i.e. the sum across 
% the rows of BLK should equal [Nu Ny].
%
% See also: balance.

% Parse Inputs
nin = nargin;
narginchk(1, 2)
[ny,nu]=size(sys);
if nin==1
    % Add single block for I/O channels
    blk = [nu ny];
end

% Get packed state-space matrix
[A,B,C,D]=ssdata(sys);
M = [A B; C D];
nx = size(A,1);

% Balance packed state-space matrix
stateblk = ones(nx,2);
allblk = [stateblk; blk];
[dMd,dr,dci] = lpvbalance(M,allblk);
if nin==1
    % Obtain similarity transformation by stripping off I/O scaling
    dr = dr/dr(end,end);
    dr = dr(1:nx,1:nx);
    
    dci = dci/dci(end,end);
    dci = dci(1:nx,1:nx);
end

% Extract balanced state matrices
Abal = dMd(1:nx,1:nx);
Bbal = dMd(1:nx,nx+1:end);
Cbal = dMd(nx+1:end,1:nx);
Dbal = dMd(nx+1:end,nx+1:end);

% Pack up data for output
% TODO PJS 5/20/2011: LPVINTERP creates an Data from PMAT state
% space matrices and then sets the various SS properties.  The method
% below (overwriting Data.a, etc) seems less prone to error since
% we don't have to worry about failing to set some properties. However,
% should we be leaving all properties untouched here?  Should StateNames
% be blown away since we are performing a state transformation?
dMd = sys;
dMd.DataPrivate.a = Abal.DataPrivate;
dMd.DataPrivate.b = Bbal.DataPrivate;
dMd.DataPrivate.c = Cbal.DataPrivate;
dMd.DataPrivate.d = Dbal.DataPrivate;


