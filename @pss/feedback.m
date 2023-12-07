function out = feedback(A,B,varargin)
% FEEDBACK  Feedback connection of two PSS objects.
%
% M = FEEDBACK(M1,M2) computes a closed-loop model M for the negative
% feedback interconnection with M1 in the forward loop and M2 in the
% feedback path, as shown in the figure below.  M is the feedback 
% interconnection at each point in the combined domains of M1 and M2. 
% To apply positive feedback, use the syntax M = FEEDBACK(M1,M2,+1).
%
%     ----->0------>[ M1 ]----------+----> 
%           |-                      |
%           |                       |
%           +-------[ M2 ]<---------- 
%
% M = FEEDBACK(M1,M2,FEEDIN,FEEDOUT,SIGN) builds a more general feedback 
% interconnection of PSS objects.  Refer to the help for 
% InputOutputModel/FEEDBACK for more details on this general interconnection 
% syntax.
%
% See also: feedback, lft, parallel, series.

out = binop(A,B,'feedback',varargin{:});

