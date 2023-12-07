%% LPVLOOPSHAPE - Parameter-varying loop-shaping synthesis
%
%  
%% Syntax
%
%    [K,GAM,INFO]=lpvloopshape(G,Wi,Wo)
%    [K,GAM,INFO]=lpvloopshape(G,Wi,Wo,Xb,Yb)
%    [K,CL,GAM,INFO]=lpvloopshape(G,...,OPT)
%
%% Description
% 
% |[K,GAM,INFO]=lpvloopshape(G,Wi,Wo)| synthesizes the parameter-varying
% controller |K|, which minimizes the induced $L_2$ norm for the shaped 
% plant |Gs=Wo*G*Wi|. |GAM| is the induced $L_2$ norm of the closed-loop system 
% from [d1,d2] to |[e1,e2].
%
%                       ^ e1                          ^ e2
%            ____       |    ____    ____     ____    |
%     +---->| Ks |-- + ---->| Wi |--| G  |-->| Wo |------->+ ----- 
%   - |      ----    ^       ----    ----     ----         ^     |
%     |              |                                     |     | 
%     |           d1 |                                  d2 |     |
%     |__________________________________________________________|
%
% |INFO| is a structure containing data from the Linear Matrix Inequalities 
% that are solved to obtain |K|. A call to |lpvloopshape| without a |basis| 
% function argument generates a controller assuming no bounds on the 
% parameter rate of variation. 
%
% |[K,GAM,INFO]=lpvloopshape(G,Wi,Wo,Xb,Yb)| computes the rate-bounded 
% controller |K|, where the rate-bounds of the independent variables of the 
% shaped pland |Gs| are incuded in the synthesis conditions. |Xb| and |Yb| are 
% |basis| objects, which describe the assumed parameter dependence of the 
% lyapunov matrices used in solving for |K|.
%
% |[K,CL,GAM,INFO]=lpvloopshape(G,...,OPT)| allows the user to pass in a 
% |lpvsynOptions| object. 
