function out = series(A,B,varargin)
% SERIES  Series connection of two UPSS objects.
%
% M = SERIES(M1,M2,OUTPUTS1,INPUTS2) connects the UPSS objects M1 and M2 in
% series. The vectors of indices OUTPUTS1 and INPUTS2 specify which 
% outputs of M1 and which inputs of M2 are connected together. M is the 
% series interconnection at each point in the combined domains of M1 and 
% M2. Refer to the help for InputOutputModel/SERIES for more details on  
% this interconnection.
%  
% If OUTPUTS1 and INPUTS2 are omitted, SERIES connects M1 and M2 in 
% cascade and returns M = M2 * M1.
%
% See also: series, append, parallel, feedback, lft.

out = binop(A,B,'series',varargin{:});

