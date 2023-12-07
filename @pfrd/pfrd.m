% PFRD   Create a parameter-varying frequency response data model
%
%   S = PFRD(Data,Domain) creates a parameter-varying frequency
%   response data model defined on an N-dimensional rectangular grid.
%   Domain is an RGRID object that  specifies the N independent variables
%   and the rectangular grid domain. Data is an N-dimensional frequency
%   response data (FRD) array.  Data(:,:,i1,...,iN) is the frequency
%   response data evaluated at the point Domain(i1,....,iN).
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 1-by-1 frequency response model defined on a 1-dimensional grid
%   IVData = linspace(2,20,10);
%   Domain = rgrid('a',IVData);
%   omeg = logspace(-1,2,30);
%   for i=1:length(IVData)
%       sys = ss(-IVData(i),IVData(i),1,0);
%       Data(1,1,i) = frd(sys,omeg);
%   end
%   S = pfrd(Data,Domain)
% 
%   % Overlay Bode plots at each independent variable
%   bode(S);
%
% See also: frd, rgrid, pgrid, pmat, pss, upmat, upss, upfrd, pstruct.

% TODO PJS 5/1/2011: This implementation requires frequency to be constant
% across the domain. Allow it to be varying?

classdef (InferiorClasses={?ss,?tf,?zpk,?frd,?pmat}) pfrd
   % Class definition for parameter-varying FRD
   
   properties (Hidden = true)
      DomainPrivate = rgrid;
      DataPrivate = frd;
   end
   
   properties (Dependent)
      Domain;
      Data;
      InputName;
      OutputName;
      Frequency;
      Parameter;
   end
   
   methods
      
      % Constructor
      function obj = pfrd(DataPrivate,Domain,varargin)
         if nargin==1
            if isa(DataPrivate,'pfrd');
               obj = DataPrivate;
            elseif isa(DataPrivate,'frd') 
               obj = pfrd(DataPrivate,rgrid);
            else
               error('Unsupported conversion to PFRD.');
            end
         elseif nargin==2 && isa(Domain,'rgrid')
            % Set properties of a pfrd
            niv = Domain.NumIV;
            szAD = size(DataPrivate);
            szad = szAD(2+niv+1:end);
            obj.DomainPrivate = rgrid(Domain,szad);
            obj.DataPrivate =  DataPrivate;
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
               error(errstr);
            end
         elseif nargin~=0
            % pfrd(Reponse,Freq,Ts) where Reponse is a PMAT
            Response = DataPrivate;
            Freq = Domain;
            if nargin==3
               Ts = varargin{3};
            else
               Ts = 0;
            end
            obj.DataPrivate = frd(Response.DataPrivate,Freq,Ts);
            obj.DomainPrivate = Response.DomainPrivate;
         end
      end
      
      % isvalid
      function [pflag,errstr] = isvalid(obj)
          % ISVALID Determine if PFRD object is valid.
         errstr = [];
         pflag = 1;
         if ~isa(obj.DomainPrivate,'rgrid')
            pflag = 0;
            errstr = ['Domain must be an rgrid object.'];
         elseif ~isa(obj.DataPrivate,'frd')
            pflag = 0;
            errstr = ['Data must be a FRD array.'];
         else
            % Data has r-by-c-by-Domain-by-ExtraDim
            szDomain = size(obj.DomainPrivate);
            szData = size(obj.DataPrivate);
            lszD = length(szDomain);
            lszS = length(szData);
            szDomain = [szDomain ones(1,lszS-2-lszD)];
            szData = [szData ones(1,lszD-lszS+2)];
            if ~( all(szDomain == szData(3:end)) || ...
                  (isempty(obj.DomainPrivate) && ndims(obj.DataPrivate)==2)   )
               pflag = 0;
               errstr = 'Dimensions of Domain and Data are incompatible.';
            end            
         end
         
         % TODO PJS 4/30/2011: Anything else to check?
         % Check for varying frequency grids?
         
         if pflag==0 && nargout==0
            error(errstr);
         end
      end
      
      % InputName
      %XXX subsref/subsasgn should work on any property of a frd
      function out = get.InputName(obj)
         out = obj.DataPrivate.InputName;
      end
      function obj = set.InputName(obj,Val)
         obj.DataPrivate.InputName = Val;
      end
      
      % OutputName
      function out = get.OutputName(obj)
         out = obj.DataPrivate.OutputName;
      end
      function obj = set.OutputName(obj,Val)
         obj.DataPrivate.OutputName = Val;
      end

      % Frequency
      function out = get.Frequency(obj)
         out = obj.DataPrivate.Frequency;
      end
      function obj = set.Frequency(obj,Val)
         obj.DataPrivate.Frequency = Val;
      end

      
      % get
      function Data = get.Data(obj)         
         DomainPrivate = obj.DomainPrivate;
         if ~isempty(DomainPrivate)
            ArrayName = DomainPrivate.ArrayName;
            idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
            % Ordering: Parameter vars, Array dims
            permidx = [find(~idx); find(idx)];
            if length(permidx) == 1
                Data = obj.DataPrivate;
            else
                Data = permute(obj.DataPrivate,permidx);
            end           
         else
            Data = obj.DataPrivate;
         end
      end
      
      % get
      function Domain = get.Domain(obj)         
         DomainPrivate = obj.DomainPrivate;
         IVNamePrivate = DomainPrivate.IVName;
         IVDataPrivate = DomainPrivate.IVData;
         IVRateBoundsPrivate = DomainPrivate.IVRateBounds;
         
         ArrayName = DomainPrivate.ArrayName;
         idx =strncmp(ArrayName,DomainPrivate.IVName,length(ArrayName));
         Domain = rgrid( IVNamePrivate(~idx), IVDataPrivate(~idx),...
                         IVRateBoundsPrivate(~idx,:));
      end
      
      % XXX should user be able to set domain parameters by hand.
      % Set Domain
      function out = set.Domain(obj,Domain)
          InData = obj.Data;
          out = pfrd(InData,Domain);
      end
      
      
      % Display
      function s=display(obj)
          % Make Header
          niv = obj.Domain.NumIV;
          S = obj.DataPrivate;
          szm = size(S);
          nfreq = length(S.Frequency);
          Ts = S.Ts;
          
          % Make header string
          [adcs,nad] = ad2char(obj.DomainPrivate);
          if nad == 0
              s1 = 'PFRD with ';
          else
              s1 = [adcs ' array of PFRDs with '];
          end
          % AH - 4/8/14 - FRDs don't keep track of the number of states.
%           s1 = [s1 int2str(ns) ' States, '];
          s1 = [s1 int2str(szm(1)) ' Outputs, '];
          s1 = [s1 int2str(szm(2)) ' Inputs, '];
          if Ts==0
              s1 = [s1 'Continuous System, '];
          else
              s1 = [s1 'Discrete System, Ts = ' num2str(Ts) ', '];
          end
          s1 = [s1 int2str(nfreq) ' Frequency points.'];
          
          % Make display string
          if niv>0
              dispStr = display(obj.Domain);
              tmp = 'The PFRD consists of the following blocks:';
              dispStr = char(tmp,dispStr(2:end,:));
          else
              dispStr = '';
          end
          s = char(s1,dispStr);
          if nargout==0
              disp(s)
              
              % TODO PJS 4/4/2011: Revisit. The code below displays the
              % matrix data if the the PMAT has only 1 IV point.
              if prod( szm(3:end) ) == 1 && ~isempty(obj)
                  %             disp([inputname(1) ' = ']);
                  %             disp(obj.Data)
                  eval([inputname(1) ' = obj.Data'])
                  
              end
          end
      end
            
%       function display(obj)
%          % Make Header
%          niv = obj.Domain.NumIV;
%          S = obj.DataPrivate;
%          szm = size(S);
%          nfreq = length(S.Frequency);
%          Ts = S.Ts;
%          if szm(1)==1
%             ioStr = ['1 Output, '];
%          else
%             ioStr = [int2str(szm(1)) ' Outputs, '];
%          end
%          if szm(2)==1
%             ioStr = [ioStr '1 Input, '];
%          else
%             ioStr = [ioStr int2str(szm(2)) ' Inputs, '];
%          end
%          if Ts==0
%             ioStr = [ioStr 'Continuous System, '];
%             %elseif Ts==-1
%             %    ioStr = [ioStr 'Discrete System, unspecified Ts, '];
%          else
%             ioStr = [ioStr 'Discrete System, Ts = ' num2str(Ts) ', '];
%          end
%          if nfreq==1
%             ioStr = [ioStr '1 Frequency point, '];
%          else
%             ioStr = [ioStr int2str(nfreq) ' Frequency points, '];
%          end
%          
%          if niv==1
%             ioStr = [ioStr '1 IV'];
%          else
%             ioStr = [ioStr int2str(niv) ' IVs'];
%          end
%          hdr = [ioStr];
%          
%          % Make header string
%          [adcs,nad] = ad2char(obj.DomainPrivate);
%          if nad==0            
%             disp(['PFRD: ' hdr]);
%          else
%             disp([adcs ' array of PFRD: ' hdr]);            
%          end       
%                            
%          % Make display string
%          if niv>0
%             dispStr = display(obj.Domain);
%          else
%             dispStr = '';
%          end
%          disp(dispStr)
%          
%          % TODO PJS 4/30/2011: Revisit. The code below displays the
%          % FRD if the the PFRD has only 1 IV point.
%          if prod( szm(3:end) ) == 1 && ~isempty(obj)
%             eval([inputname(1) ' = S'])
%          end
%       end
      
   end % end of methods
   
   methods (Access = private)
      % privatesize
      function varargout = privatesize(s,arg2)
         out = size(s.DataPrivate);
         if length(out)==3
            out = [out 1];
         end

         if nargin==1
            arg2 = nan;
         end

         nout = nargout;
         if nout==2
             % Direct call to csize with nout = 2 put sthe product of the column
             % and IV dims into the second output. Instead SIZE with two 
             % output args should return only the row/col dimensions.
             tmpout = csize(out,arg2,nargin,3);
             varargout = {tmpout{1}, tmpout{2}};
         else
             varargout = csize(out,arg2,nargin,nout);
         end
      end
  end
   
end % end of classdef





