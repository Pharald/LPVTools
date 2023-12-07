%% UPMAT - Uncertain parameter-varying matrix
%
%  
%% Syntax
%
% |M = upmat(Data,Domain)|
%
%% Description
%
% |M = upmat(Data,Domain)| creates an uncertain parameter-varying matrix
% defined on an N-dimensional rectangular grid. |Domain| is an |rgrid| object
% that specifies the N independent variables and the rectangular grid
% domain.  |Data| is an (N+2) dimensional double array that specifies the
% matrix data. |Data(:,:,i1,...,iN)| is the value of the matrix evaluated
% at the point |Domain(i1,....,iN)|.
% 
%% Example

% Create a 2-by-2 matrix defined on a 1-dimensional grid
IVData = linspace(-2,2,20);
Domain = rgrid('a',IVData);
au = ureal('au',-2.3);
for i=1:length(IVData)
    Data(1:2,1:2,i) = [1 au*IVData(i); IVData(i)^2 cos(IVData(i))];
end
M = upmat(Data,Domain)