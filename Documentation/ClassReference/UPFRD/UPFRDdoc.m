%% UPFRD - Uncertain parameter-varying frequency response data model
%
%  
%% Syntax
%
% |S = upfrd(Data,Domain)|
%
%% Description
%
% |S = upfrd(Data,Domain)| creates an uncertain parameter-varying
% frequency response data model defined on an N-dimensional rectangular
% grid. |Domain| is an |rgrid| object that  specifies the N independent
% variables and the rectangular grid domain. |Data| is an N-dimensional
% uncertain frequency response data (|ufrd|) array.
% |Data(:,:,i1,...,iN)| is the uncertain frequency response data evaluated
% at the point |Domain(i1,....,iN)|.
% 
%% Example


% Create a 1-by-1 UFRD defined on a 1-dimensional grid
IVData = linspace(2,20,10);
Domain = rgrid('a',IVData);
omeg = logspace(-1,2,30);
unc = ureal('unc',10) 
usys = rss(1,1,2,10)*unc;
Data = ufrd(usys,omeg);
S = upfrd(Data,Domain)

%%
% Overlay Bode plots at each independent variable
bode(S);