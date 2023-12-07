function han = lpvplot(varargin)
% LPVPLOT  Plot PMAT data as a function of the independent variable.
%
% LPVPLOT plots the data contained in a PMAT against the independent 
% variable it depends on. The syntax is identical to the UPLOT command 
% in the robust control toolbox, except the data are PMATs.
%
% LPVPLOT(PLOT_TYPE,PMAT1,PMAT2,PMAT3, ...)  plots the value of PMAT1,
% PMAT2, PMAT3, ... against the independent variable. PMAT1, PMAT2, 
% PMAT3, ... can only depend on a single independent variable, and they
% must all depend on the same independent variable.
%
% The (optional) PLOT_TYPE argument must be one of:
%
%   'iv,d'       matin .vs. independent variable (default option)
%   'iv,m'       magnitude .vs. independent variable
%   'iv,lm'      log(magnitude) .vs. independent variable
%   'iv,p'       phase .vs. independent variable
%   'liv,d'      matin .vs. log(independent variable)
%   'liv,m'      magnitude .vs. log(independent variable)
%   'liv,lm'     log(magnitude) .vs. log(independent variable)
%   'liv,p'      phase .vs. log(independent variable)
%   'nyq'        real .vs. imaginary  (parametrized by indep variable)
%   'nic'        Nichols chart
%   'bode'       Bode magnitude and phase plots
%
% LPVPLOT(PLOT_TYPE,PMAT1,'linetype1',PMAT2,'linetype2',...) enables the
% user to set the linetype for the data associated with each PMAT input.
%
%See also: BODE, LOGLOG, PLOT, NICHOLS, NYQUIST, SEMILOGX, SEMILOGY, SIGMA.

% TODO PJS 4/4/2011: Revisit

% Find PMATs with nonempty ranges and do error checking
nin = nargin;
ivdrange = [];
ivassigned = false;
for i=1:nin
   v = varargin{i};
   if isa(v,'pmat')
      niv = v.Domain.NumIV;
      ivn = v.Domain.IVName;  
      if ivassigned == false
        ivnfirst = ivn;
        ivassigned = true;
      end
      ivd = v.Domain.IVData;         
      if isempty(ivdrange)
         if niv>1
            error('All PMATs must have single independent variable.');
         elseif niv==1 && ~isempty(ivn{1})
            ivrange = [ivd{1}(1) ivd{1}(end)];
         end
      else
         if ~isequal(ivnfirst,ivn) 
            error('All PMATs must have same independent variables');
         end
         ivrange(1) = min(ivrange(1), ivd{1}(1));
         ivrange(2) = max(ivrange(2), ivd{1}(end));
      end      
   end
end

% Replace PMATs with double data
newin = cell(0);
for i=1:nin
   v = varargin{i};
   if isa(v,'pmat')
      niv = v.Domain.NumIV;
      ivd = v.Domain.IVData;         
      Data = v.Data;
      szm = [privatesize(v) 1];
      npts = prod(szm(3:end));
      if niv==0
         ivd = {ivrange};
         Data = repmat(Data,[1 1 2]);
         npts = 2;
      end      
      Data = reshape(Data,[szm(1)*szm(2) npts]).';
      newin = [newin {ivd{1}, Data}];
   else
      newin = [newin varargin{i}];
   end
end
if nargout==0
   uplot(newin{:});
   xlabel(ivnfirst);
else
   han = uplot(newin{:});
   xlabel(ivnfirst);
end
