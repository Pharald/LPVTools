function out = loopsens(varargin)
% LOOPSENS  Feedback loop sensitivity functions for PLFT objects
%
% LOOPTRANSFER = LOOPSENS(P,C,DOMAIN) constructs the loop transfer 
% functions for the multivariable feedback loop consisting of the PLFT C in 
% negative feedback with the PLFT P, as shown in the diagram below.
%
%     ----->0------>[ P ]----------+---->
%           |-                     |
%           |                      |
%     <-----+-------[ C ]<---------0<----
%
% LOOPTRANSFER is a PSTRUCT containing fields for all input/output 
% transfer functions of the feedback loop at each point in DOMAIN 
% It also has fields Poles and Stable that contain the closed-loop 
% poles and a flag for the stability of the feedback loop. Additional 
% details on the fields of LOOPTRANSFER can be found in the documentation 
% for DynamicSystem/LOOPSENS.
%
% See also: loopsens, robuststab, robustperf, wcgain.

DOM = varargin{end};
varargin(end) = [];

if ~isa(DOM,'rgrid')
    error(['Last input argument must be an RGRID object that specifies'...
           ' the parameter domain grid.'])
end

for i = 1:numel(varargin)
    if isa(varargin{i},'plft')
        varargin{i} = lft2grid(varargin{i},DOM);
    end
end

out = loopsens(varargin{:});
end