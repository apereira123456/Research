x=[0:0.01:10]';

Spectrum = normpdf(x,3,0.5)+5*normpdf(x,7,1)+0.2*x+0.01*x.^2;

[Base, Corrected_Spectrum]=baseline(Spectrum);

plot(x,Spectrum,'b-',x,Base,'r--',x,Corrected_Spectrum,'g-.');

legend('Initial Spectrum','Baseline','Corrected Spectrum');