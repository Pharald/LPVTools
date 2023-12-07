function [A,B] = switchyard(A,B)
% SWITCHYARD   Utility function for LPV object interaction. 


% Handle: double, ss, frd, tf, zpk, umat,uss,ufrd, pmat,pfrd,pss,
%         upmat, upss, upfrd
objcell = {'upfrd', 'upss', 'upmat', ...  % 1,2,3
           'pfrd', 'pss', 'pmat', ... % 4,5,6
           'frd', 'ss', 'tf', 'zpk', 'double', ... % 7,8,9,10,11
           'ufrd', 'uss', 'umat'}; % 12,13,14
fcell = {'upfrd', 'ufrd', 'frd', 'pfrd'};
ucell = {'upfrd', 'upss', 'upmat', 'ufrd', 'uss', 'umat'};
pcell = {'upfrd', 'upss', 'upmat', 'pfrd', 'pss', 'pmat'};
scell = {'upss', 'pss', 'ss', 'tf','zpk','uss'};

% Lift pgrid object immediately to a pmat
if isa(A,'pgrid');
    A = pmat(A);
end
if isa(B,'pgrid');
    B = pmat(B);
end

% Is one a FRD and the other not?
Aftype = ismember(class(A),fcell);
Bftype = ismember(class(B),fcell);

% Promote uncertain atoms to objects
A = LOCALLiftAtom(A);
B = LOCALLiftAtom(B);

% TODO: Frequency Units
if xor(Aftype,Bftype) % One is FRD and other is not
   % Switch order so that A is the FRD   
   flg = (Aftype == 0 && Bftype ==1);
   [A,B] = LOCALSwitchArgs(A,B,flg);
   
   % Lift B to be FRD
   Boloc = find(strcmp(class(B),objcell));
   omega = A.Frequency;
   nomega = length(omega);
   switch Boloc
      case {2,3} % up
         B = upfrd(B,omega);
      case {5,6} % p
         B = frd( pss(B), omega );
      case {8,9,10,11} % normal
         B = frd(pss(B),omega);
      case {13,14} % u
         B = ufrd(upss(B),omega);
   end
   [A,B] = LOCALSwitchArgs(A,B,flg);
end

% Is one a Uncertain and the other not?
Autype = ismember(class(A),ucell);
Butype = ismember(class(B),ucell);

if xor(Autype,Butype) % One is Uncertain and other is not
   % Switch order so that A is the Uncertain   
   flg = (Autype == 0 && Butype ==1);
   [A,B] = LOCALSwitchArgs(A,B,flg);
   
   % Lift B to be Uncertain
   Boloc = find(strcmp(class(B),objcell));
   switch Boloc
      case {4}
         B = upfrd(B);
      case {5}
         B = upss(B);
      case {6}
         B = upmat(B);
      case {7}
         B = ufrd(B);
      case {8,9,10}
         B = uss(B);
      case {11}
         B = umat(B);
   end
   [A,B] = LOCALSwitchArgs(A,B,flg);
end


% Is one a Parameter-Varying and the other not?
% Note: At least one of A and b must be parameter-varying to get here.
Aptype = ismember(class(A),pcell);
Bptype = ismember(class(B),pcell);
  
if xor(Aptype,Bptype) % One is Parameter-Varying and other is not
   % Switch order so that A is the Parameter-Varying   
   flg = (Aptype == 0 && Bptype ==1);
   [A,B] = LOCALSwitchArgs(A,B,flg);

   % Lift B to be Parameter-Varying
   Boloc = find(strcmp(class(B),objcell));
   switch Boloc
      case {7}
         B = pfrd(B);
      case {8,9,10}
         B = pss(B);
      case {11}
         B = pmat(B);
      case {12}
         B = upfrd(B);
      case {13}
         B = upss(B);
      case {14}
         B = upmat(B);
   end
   [A,B] = LOCALSwitchArgs(A,B,flg);
end


% Is one an S-type and the other not?
Astype = ismember(class(A),scell);
Bstype = ismember(class(B),scell);
if xor(Astype,Bstype) % One is an S-type and other is not
   % Switch order so that A is the S-type
   flg = (Astype == 0 && Bstype ==1);
   [A,B] = LOCALSwitchArgs(A,B,flg);
      
   % Lift B to be S-type
   % Note: At this point A and B are both parameter-varying (because
   % the switchyard only gets called if A or B is parameter-vyaring)
   % Thus cases 11 and 14 currently cannot occur.
   Boloc = find(strcmp(class(B),objcell));
   switch Boloc
      case {3} % up
         B = upss(B);
      case {6} % p
         B = pss(B);
      case {11} % normal
         B = ss(B);
      case {14} % u
         B = uss(B);
   end
   
   [A,B] = LOCALSwitchArgs(A,B,flg);   
end

% Local function to lift atoms to uncertain objects
function A = LOCALLiftAtom(A)

if isuncertain(A)
    if isa(A,'ucomplex') || isa(A,'ucomplexm') || isa(A,'ureal');
        A = umat(A);
    elseif isa(A,'udyn') || isa(A,'ultidyn')
        A = uss(A);
    end
end


% Local function to switch arguments
function [A,B] = LOCALSwitchArgs(A,B,flg)

if flg
   tmp = A;
   A = B;
   B = tmp;
end


