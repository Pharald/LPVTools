function out = parallel(A,B,varargin)
% PARALLEL  Parallel connection of two PSS objects.
%
% M = PARALLEL(M1,M2,IN1,IN2,OUT1,OUT2) connects the PSS objects M1 and M2
% in parallel. The inputs specified by IN1 and IN2 are connected and the
% outputs specified by OUT1 and OUT2 are summed. M is the parallel 
% interconnection at each point in the combined domains of M1 and M2.
% Refer to the help for InputOutputModel/PARALLEL for more details on this 
% interconnection.
%
% If IN1,IN2,OUT1,OUT2 are omitted, PARALLEL forms the standard parallel
% interconnection of M1 and M2 and returns M = M1 + M2.
%
% See also: parallel, append, series, feedback, lft.

out = binop(A,B,'parallel',varargin{:});

