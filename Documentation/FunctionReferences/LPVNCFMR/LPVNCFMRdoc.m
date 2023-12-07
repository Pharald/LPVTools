%% LPVNCFMR - Contractive coprime factor model reduction of a |pss|
%  
%% Syntax
%
%    [Pred,INFO] = lpvncfmr(P)
%    [Pred,INFO] = lpvncfmr(P,ORDER)
%
%% Description
%
% |lpvncfmr| performs balanced truncation model reduction through contractive
% coprime factorization of a |pss|.
% 
% |[Pred,INFO] = lpvncfmr(P,ORDER)| finds a balanced contractive
% coprime factorization of the LPV system |P| (analogous to a normalized
% coprime factorization for LTI systems), and performs a balanced 
% truncation to remove those states that contribute the least to the 
% input-output mapping of the balanced LPV system. If |P| has M states, then 
% the reduced order system Pred will have |ORDER| states, with M - |ORDER| states 
% removed using the balanced truncation. |INFO.hsv| contains a vector of 
% singular values describing the input-output mapping of |SYSB| (comparable 
% to Hankel singular values for LTI systems).
%
% |[Pred,INFO] = lpvncfmr(P)| computes the balanced contractive coprime 
% factorization of the LPV system |P|, and outputs it as |Pred|. 
% This is equivalent to the call |[Pred,INFO] = LPVNCFMR(P,Nx)| for a 
% system |P| with |Nx| states.