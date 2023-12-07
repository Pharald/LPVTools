%% PFRD - Parameter-varying frequency response data model
%
%  
%% Syntax
%
% |S = pfrd(Data,Domain)|
%
%% Description
%
% |S = pfrd(Data,Domain)| creates a parameter-varying frequency
% response data model defined on an N-dimensional rectangular grid.
% |Domain| is an |rgrid| object that  specifies the N independent variables
% and the rectangular grid domain. |Data| is an N-dimensional frequency
% response data (|frd|) array.  |Data(:,:,i1,...,iN)| is the frequency
% response data evaluated at the point |Domain(i1,....,iN)|.
% 
%% Example

% Create a 1-by-1 frequency response model defined on a 1-dimensional grid
IVData = linspace(2,20,10);
Domain = rgrid('a',IVData);
omeg = logspace(-1,2,30);
for i=1:length(IVData)
  sys = ss(-IVData(i),IVData(i),1,0);
  Data(1,1,i) = frd(sys,omeg);
end
S = pfrd(Data,Domain) 

%%
% Overlay Bode plots at each independent variable
bode(S);