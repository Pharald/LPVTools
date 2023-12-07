%% LPVBALANCE - Diagonal scaling of a |pmat| or |pss| 
%
%  
%% Syntax
%
%    [B,D] = lpvbalance(A)
%    [B,DR,DC] = lpvbalance(A,BLK)
%
%% Description
% 
% |lpvbalance| computes a diagonal scaling for LPV objects to improve their
% numerical conditioning. The algorithm used to accomplish this uses a 
% generalized version of Osborne's iteration.
%
% *|lpvbalance| for |pmat|*
%
% |[B,D] = lpvbalance(A)| computes a single diagonal similarity transformation
% to improve the conditioning of the N-by-N |pmat| |A| at all points in the
% domain.  The transformation |D| is computed using a generalized version of
% Osborne's iteration. |D| is returned as an N-by-N diagonal, |double| matrix.
% The scaled |pmat| |B| is |D*A*inv(D)|.  This applies the transformation |D| at
% each point in the domain of |A|.
%
% |[B,DR,DC] = lpvbalance(A,BLK)| computes structured, diagonal transformations
% for the N-by-M |pmat| |A|. The scaled |pmat| |B| is |DR*A*DC| where |DR| is an
% N-by-N matrix and and |DC| is an M-by-M matrix. |BLK| is a K-by-2 matrix
% that specifies the block partitioning dimensions of |DR| and |DC|. If
% |BLK = [c1 r1; ... ; ck rk]| then |DR| and |DC| are partitioned as:
%
%    DR = blkdiag( d1*I_r1, ...,  dk*I_rk )
%    DC = blkdiag( (1/d1)*I_c1, ...,  (1/dk)*I_ck )
%
% where the notation |I_r1= eye(r1)|, represents the r1-by-r1 identity matrix.
% The block partitioning must be consistent with the dimensions of |A|, i.e.
% the sum across the rows of |BLK| should equal [M N].
% 
% 
% *|lpvbalance| for |pss|*
%
% |[P,D] = lpvbalance(S)| computes a single diagonal similarity transformation 
% to improve the conditioning of the Ny-by-Nu |pss| |S| at all points in the
% domain. The transformation |D| is computed using a generalized version of
% Osborne's iteration. |D| is returned as an Nx-by-Nx diagonal, double 
% matrix where Nx is the state dimension of |S|. The scaled |pss| |P| is 
% obtained by applying the similarity transformation |D| at each point in 
% the domain of |S|.
%
% |[P,DR,DC] = lpvbalance(S,BLK)| computes diagonal transformations applied to
% the states and input/output channels. The state matrices of the scaled 
% |pss| |P| are obtained from |DR*[A B; C D]*DC| where |A|, |B|, |C|, |D|
% are the state matrices of |S|. |BLK| is a K-by-2 matrix that specifies 
% the partitioning dimensions of |DR| and |DC|. If |BLK = [c1 r1; ... ; ck rk]| 
% then |DR| and |DC| are partitioned as:
%
%    DR = blkdiag( D, d1*I_r1, ...,  dk*I_rk )
%    DC = blkdiag( inv(D), (1/d1)*I_c1, ...,  (1/dk)*I_ck )
%
% where |D| is a diagonal, Nx-by-Nx matrix, and the notation |I_r1= eye(r1)|, 
% represents the r1-by-r1 identity matrix. The block partitioning must be 
% consistent with the input/output dimensions of |S|, i.e. the sum across 
% the rows of |BLK| should equal [Nu Ny].