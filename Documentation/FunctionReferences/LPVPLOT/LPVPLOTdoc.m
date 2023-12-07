%% LPVPLOT - Plot PMAT data as a function of the independent variable.
%
%  
%% Syntax
%
%    lpvplot(PMAT)
%    lpvplot(PLOT_TYPE,PMAT1,PMAT2,PMAT3, ...)
%    lpvplot(PLOT_TYPE,PMAT1,'linetype1',PMAT2,'linetype2',...)
%
%% Description
% 
%
% |lpvplot| plots the data contained in a |pmat| against the independent 
% variable it depends on. The syntax is identical to the |uplot| command 
% in the Robust Control Toolbox, except the data are |pmat|.
%
% |lpvplot(PLOT_TYPE,PMAT1,PMAT2,PMAT3, ..., PMATN)|  plots the value of |PMAT1|,
% |PMAT2|, |PMAT3|, ..., |PMATN| against the independent variable. 
% |PMAT1|, |PMAT2|, |PMAT3|, ... can only depend on a single independent 
% variable, and they must all depend on the same independent variable.
%
% The (optional) PLOT_TYPE argument must be one of:
%
% * |'iv,d'|       matin .vs. independent variable (default option)
% * |'iv,m'|       magnitude .vs. independent variable
% * |'iv,lm'|      log(magnitude) .vs. independent variable
% * |'iv,p'|       phase .vs. independent variable
% * |'liv,d'|      matin .vs. log(independent variable)
% * |'liv,m'|      magnitude .vs. log(independent variable)
% * |'liv,lm'|     log(magnitude) .vs. log(independent variable)
% * |'liv,p'|      phase .vs. log(independent variable)
% * |'nyq'|        real .vs. imaginary  (parametrized by indep variable)
% * |'nic'|        Nichols chart
% * |'bode'|       Bode magnitude and phase plots
%
% |lpvplot(PLOT_TYPE,PMAT1,'linetype1',PMAT2,'linetype2',...)| enables the
% user to set the linetype for the data associated with each |pmat| input.




