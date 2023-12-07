% function [N,W]=covertd(taud,ordp,ordw)
%
% Cover a time delay exp(-j*om*taud) with an uncertainty set
%      U := { N(s)+Delta(s) :  |Delta(jw)| <= |W(jw)| }
% The nominal approximation N(s) is a Pade approximation of
% order ordp.  The weight W(s) overbounds the approximation error
% at each frequency and is of order ordw.
function [N,W,Efr]=covertd(taud,ordp,ordw)

% Construct Delay
D = ss(1);
D.InputDelay = taud;

% Nominal Pade Model and Approximation Error
N = pade(D,ordp);
E = D-N;

% Response of error on frequency grid
% XXX Freq grid can be an input option
Nw = 1e3;
om = logspace(-4,2,Nw)/taud;
Efr = freqresp(E,om);

% Smooth error to improve numerical conditioning of FITMAGFRD
% XXX The high and low frequency smoothing parameters can be options.
hfgain = 2;
hftol = 0.05;
idx = min(find(abs(Efr)>hfgain-hftol));
Efr(idx:end) = hfgain;

lfgain =  1e-4;
idx = max(find(abs(Efr)<lfgain));
Efr(1:idx) = lfgain;

% Weight Overbound
Efr = frd(Efr,om);
W =  fitmagfrd(Efr,ordw,0,[],1);

