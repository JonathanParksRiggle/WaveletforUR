% function [fs,tau,qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice)
% 
% Calculates scale and period arrays for wavelet transform
% Script written by Tanya Leise, Department of Mathematics, Amherst College
function [fs,tau,qscaleArray] = CalcScaleForCWT(shortestperiod,longestperiod,T,NstepsPerHr,nvoice)

scale=floor(2*pi*T/longestperiod);
noctave = ceil(log2(2*pi*T/scale/shortestperiod));
nscale  = nvoice .* noctave;
kscale  = 1; 

qscaleArray=zeros(nscale,1);
for jo = 1:noctave,
    for jv = 1:nvoice,
        qscale = scale .* (2^(jv/nvoice));
        qscaleArray(nscale-kscale+1)=qscale;
        kscale = kscale+1;
    end
    scale  = scale .*2;
end

tau=(2*pi)*T./qscaleArray;
fs=2*pi./(tau*NstepsPerHr);