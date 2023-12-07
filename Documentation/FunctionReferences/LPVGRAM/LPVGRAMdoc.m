%% LPVGRAM - Compute Gramians for |pss| objects
%
%  
%% Syntax
%
%    Wc = lpvgram(SYS,'c')
%    Wo = lpvgram(SYS,'o')
%    W = lpvgram(SYS,OPTION,WEIGHT)
%    W = lpvgram(SYS,...,INVERT)
%
%% Description
% 
% |Wc = lpvgram(SYS,'c')| computes the controllability gramian of the |pss| |sys|.
% The output |Wc| is a constant |double| matrix, which satisfies the 
% LMI: |A*Wc+Wc*A' +B*B' < 0| at each point in the domain of |SYS|, where |A| is 
% the state matrix of |SYS| and |B| is its input matrix.
%
% |Wo = lpvgram(SYS,'o')| computes the observability gramian of the |pss| |SYS|. 
% The output |Wo| is a constant |double| matrix, which satisfies the 
% LMI: |A'*Wo+Wo*A +C'*C < 0| at each point in the domain of |SYS|, where |A| is 
% the state matrix of |SYS| and |C| is its output matrix.
%
% |W = lpvgram(SYS,OPTION,WEIGHT)| applies a matrix weighting |WEIGHT| when 
% solving for the gramian. For a controllabilty gramian the LMI becomes: 
% 
%                 A*WEIGHT*Wc+Wc*WEIGHT*A' +B*B' < 0
% 
% For a observability gramian the LMI becomes: 
% 
%                 A'*WEIGHT*Wo+Wo*WEIGHT*A +C'*C < 0
% 
% If no |WEIGHT| is specified a default value of |eye(size(A))| is used, and
% the resulting gramian is diagonal.
%
% |W = lpvgram(SYS,...,INVERT)| provides an alternative implementation of the
% algorithm which solves for the gramians. If |INVERT| is |True| the LMI
% conditions are changed to solve for the inverse of the gramians, which
% can improve the accuracy of the solution for certain systems. The default 
% implementation assumes |INVERT=FALSE|.