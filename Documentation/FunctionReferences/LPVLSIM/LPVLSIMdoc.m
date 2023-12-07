%% LPVLSIM - Simulate the time response of a LPV system
%
%  
%% Syntax
%
%    [Y,T,X,U,TRAJ] = lpvlsim(G,PTRAJ,UIN,TIN)
%    [Y,T,X,U,TRAJ] = lpvlsim(G,PTRAJ,UIN,TIN,X0)
%
%% Description
%
% |[Y,T,X,U,TRAJ] = lpvlsim(G,PTRAJ,UIN,TIN)| simulates the time-response of 
% the system |G|, subject to the input signal defined by |UIN| and |TIN|, and the 
% parameter tracjetory defined in |PTRAJ|. |G| is a |pss| with |Ny| outputs, |Nx| states,
% |Nu| inputs, and |N| independent variables |IVName1,...,IVNameN|. |TIN| is a sorted 
% column vector of time values, and |UIN| is a |length(TIN)-by-Nu| matrix of 
% corresponding inputs. |PTRAJ| is a struct which defines the time-variation 
% of the parameters (independent variables) in |G|. The field |PTRAJ.time| 
% contains a sorted row vector of time-values. |PTRAJ| must also have a field 
% for each independend variable in |G|, such that |PTRAJ.IVName1, ... ,PTRAJ.IVName| 
% each contain a row vector of parameter trajectories corresponding to 
% |PTRAJ.time|. |Y| is a |length(T)-by-NY| matrix whose columns correspond to
% the outputs of |G|, |X| is a |length(T)-by-Nx| matrix whose columns 
% correspond to the state trajectories of |G|, |U| is a |length(T)-by-Nu| matrix 
% whose columns correspond to the inputs of |G|, and |T| is a column vector of 
% time values corresponding to |Y|, |X| and |U|. |TRAJ| contains the corresponding 
% parameter trajectories. 
%
% |[Y,T,X,U,TRAJ] = lpvlsim(G,PTRAJ,UIN,TIN,X0)| simulates the time-response 
% of the system |G| starting from the initial condition |X0|.

