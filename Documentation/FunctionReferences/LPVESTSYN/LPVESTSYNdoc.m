%% LPVESTSYN - Synthesize a parameter-varying estimator for LPV systems
%
%  
%% Syntax
%
%    [L,GAM,INFO] = lpvestsyn(P)
%    [L,GAM,INFO] = lpvestsyn(P,'L2')
%    [L,GAM,INFO] = lpvestsyn(P,'LQG')
%    [L,GAM,INFO] = lpvestsyn(P,Yb)
%    [L,GAM,INFO] = lpvestsyn(P,Yb,'L2')   
%    [L,GAM,INFO] = lpvestsyn(P,Yb,'LQG') 
% 
%% Description
% 
% |[L,GAM,INFO] = lpvestsyn(P,'L2')| computes a parameter-varying state  
% estimator for the parameter-varying system |P|. |L| takes in all the outputs 
% of |P| and outputs an estimate of the states of |P|. |L| is the constant   
% estimation matrix for the plant |P|, which minimizes the induced $L_2$ norm 
% of the error  in the state-estimate. |GAM| is the minimum $L_2$ norm achived 
% by |L|. |INFO| is a struct with additional data. 
%
% |[L,GAM,INFO] = lpvestsyn(P,'LQG')| computes the constant estimation matrix 
% for the plant |P|, which minimizes the stochastic LPV bound on the state
% estimation error. The stochastic LPV bound on the state estimation error  
% is defined as the expected value of the average instantaneous power of 
% error signal, assuming system inputs are zero mean, white-noise processes  
% with unit intensity.
%
% |[L,GAM,INFO] = lpvestsyn(P,Yb,ALG)| performs a rate-bounded synthesis. 
% |Yb| is a |basis| object specifying the basis functions to be used in the 
% synthesis. |ALG| can be either |'L2'| or |'LQG'|. A call without the |ALG|
% argument is equivalent to |[L,GAM,INFO] = lpvestsyn(P,Yb,'L2')|
%