% B4C
[FitResults,FitError] = peakfit(signal,0,0,8,1,0,1,app.first_guess,0);

% SiC
[~,peak1_start] = min(abs(signal(:,1) - 705));
[~,peak1_end] = min(abs(signal(:,1) - 875));

[~,peak2_start] = min(abs(signal(:,1) - 875));
[~,peak2_end] = min(abs(signal(:,1) - 1053));

figure
[FitResults1,FitError1,Baseline1] = peakfit(signal(peak1_start:peak1_end,:),0,0,4,1,0,1,app.first_guess.sic1,1);

figure
[FitResults2,FitError2] = peakfit(signal(peak2_start:peak2_end,:),0,0,1,1,0,1,app.first_guess.sic2,1);

function err = fitBWF(lambda,t,y,shapeconstant)
%   Fitting function for multiple Breit-Wigner-Fano.
% T. C. O'Haver (toh@umd.edu),  2014.
global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
A = zeros(length(t),round(length(lambda)/2));
for j = 1:length(lambda)/2,
    A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
end
if AUTOZERO==3,A=[ones(size(y))' A];end
if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
z = A*PEAKHEIGHTS;
if LOGPLOT,
    err = norm(log10(z)-log10(y)');
else
    err = norm(z-y');
end
% ----------------------------------------------------------------------
function g = BWF(x,pos,wid,m)
% BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
% pos=position; wid=width; m=Fano factor
%  T. C. O'Haver, 2014
y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
% y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
g=y./max(y);
end