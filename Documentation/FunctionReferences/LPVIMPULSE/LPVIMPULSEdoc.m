%% LPVIMPULSE - Parameter dependent impulse response for LPV systems
%
%  
%% Syntax
%
%    lpvimpulse(SYS,PTRAJ)
%    [Y,T,X,U,PTRAJOUT] = LPVIMPULSE(SYS,PTRAJ)
%    [Y,T,X,U,PTRAJOUT] = lpvimpulse(SYS,ptraj,TFINAL)
%    [Y,T,X,U,PTRAJOUT] = lpvimpulse(SYS,ptraj,T)
% 
%
%% Description 
% 
% |[Y,T,X,U,PTRAJOUT] = lpvimpulse(SYS,PTRAJ)| computes the parameter dependent
% impulse response of |SYS|. |SYS| is a LPV system with |Ny| outputs, |Nx| states,
% |Nu| inputs, and |N| independent variables |IVName1,...,IVNameN|. |PTRAJ| is a 
% struct which defines the time-variation of the parameters (independent 
% variables) in |SYS|. The field |PTRAJ.time| contains a sorted row vector of 
% time-values. |PTRAJ| must also have a field for each independend variable 
% in |SYS|, such that |PTRAJ.IVName1, ... ,PTRAJ.IVNameN| each contain a row 
% vector of parameter trajectories corresponding to |PTRAJ.time|. 
% The output |Y| is a |length(T)-by-NY-by-Nu| matrix such that |Y(:,i,j)|  
% corresponds to the i-th output of |SYS| due to a step command in the j-th 
% input channel. Similarly |X| is a |length(T)-by-Nx-by-Nu| matrix describing
% the state trajectories of |SYS|, |U| is a |length(T)-by-Nu-by-Nu| matrix 
% describing the trajectory of the inputs to |SYS|, and |T| is a column vector 
% of time values corresponding to |Y|, |X| and |U|. |PTRAJOUT| contains the 
% corresponding parameter trajectories.
%
% |lpvimpulse(SYS,PTRAJ)| generates plots of the parameter dependent impulse 
% response of |SYS|.
% 
% |[Y,T,X,U,PTRAJOUT] = lpvimpulse(SYS,ptraj,TFINAL)| simulates the impulse
% response up to the time |TFINAL|.
%
% |[Y,T,X,U,PTRAJOUT] = lpvimpulse(SYS,ptraj,T)| simulates the impulse response
% using a user supplied time vector |T|.
