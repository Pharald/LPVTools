function out = idf(ivname,values,ratebound)
%IDF    Creates a function whose value is the argument.
%   F = IDF(VALUES) creates a PMAT 1x1 matrix, with
%   matrix value equal to the independent variable (ie.,
%   the identity function).  The values of the independent
%   variable are VALUES.  The independent variable's name
%   is automatically set to the name of the variable VALUES
%   from the caller's workspace.
%
%   F = IDF(IVNAME,VALUES) is the same, but the name
%   of the independent variable is IVNAME.
%
%   F = IDF(IVNAME,VALUES,RATEBOUNDS) specifies the maximum and
%   minimum rate bounds using the 1-by-2 vector RATEBOUNDS.
%
%   IDF(VALUES) is the same, but the result overwrites the
%   variable VALUES in the caller's workspace.
%
%   IDF(IVNAME,VALUES,RATEBOUNDS) is the same, but the result is
%   assigned in the caller's workspace with variable name IVNAME.
%
% Example:
%   t = linspace(0,10,300);
%   idf(t)
%   lpvplot(sqrt(t)*sin(t))

% NOTE PJS: Revisit

% Check # of input/output arguments
nin = nargin;
error(nargchk(1, 3, nin, 'struct'))
nout = nargout;
error(nargchk(0, 1, nout, 'struct'))

if nin==1
    values =ivname;
    ivname = inputname(1);
    ratebound = [-inf inf];
elseif nin == 2
    ratebound = [-inf inf];
end

if ~( isa(values,'double') && isvector(values) )
    error('VALUES must be a double vector');
elseif ~( isa(ivname,'char') && size(ivname,1)==1 && ndims(ivname)==2 )
    error('IVNAME must be a character string');
elseif ~(  isa(ratebound,'double') && numel(ratebound)==2 )
    error('RATEBOUND must be a 1-by-2 double');
else
    ratebound = sort(ratebound(:))';
    values = sort(values(:));
    lv = length(values);
    domain = rgrid(ivname,values,ratebound);
    tmp = pmat(reshape(values,[1 1 lv]),domain);
    if nargout==0
        assignin('caller',ivname,tmp);
        evalin('caller',['display(' ivname ')']);
    else
        out = tmp;
    end
end
