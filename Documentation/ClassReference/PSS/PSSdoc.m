%% PSS - Parameter-varying state-space model
%
%  
%% Syntax
%
% |S = pss(Data,Domain)| 
%
%% Description
%
% |S = pss(Data,Domain)| creates a parameter-varying state-space
% model defined on an N-dimensional rectangular grid. |Domain| is an |rgrid|
% object that  specifies the N independent variables and the rectangular
% grid domain. |Data| is an N-dimensional state-space array that
% specifies the state space data. |Data(:,:,i1,...,iN)| is the model
% evaluated at the point |Domain(i1,....,iN)|.
% 
%% Example
% 

% Create a 1-by-1 state-space model defined on a 1-dimensional grid
IVData = linspace(2,20,10);
Domain = rgrid('a',IVData);
for i=1:length(IVData)
  Data(1,1,i) = ss(-IVData(i),IVData(i),1,0);
end
S = pss(Data,Domain)

%%
% Overlay Bode plots at each independent variable
bode(S)
