%% PSTRUCT - Parameter-varying structure
%
%  
%% Syntax
%
% |M = pstruct(Data,Domain)|
%
%% Description
%
% |M = pstruct(Data,Domain)| creates a parameter-varying structure defined on
% an N-dimensional rectangular grid. |Domain| is an |rgrid| object that
% specifies the N independent variables and the rectangular grid domain.
% |Data| is an (N+2) dimensional structured array that specifies the
% data. |Data(:,:,i1,...,iN)| is the value of the struct array evaluated
% at the point |Domain(i1,....,iN)|.
% 
% Note: Use |M.FieldName| to access the field named 'FieldName' in M.
% If possible, the content of the field is returned as an
% object (e.g. pmat, pss), otherwise it is returned as a cell array.
% 

