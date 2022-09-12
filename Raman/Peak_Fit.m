signal = [Wavenumber, Intensity];

% [FitResults,FitError] = peakfit(signal,center,window,NumPeaks,peakshape,extra, NumTrials, start, BaselineMode, fixedparameters, plots);
% center and window are set to 0 to fit the entire signal
% NumPeaks is 6 for B4C
% peakshape is 1 for unconstrained Gaussian
% extra might be 0 or 1
% NumTrials is 1
% start specifies the firstguess in the form [position1 width1 position2 width2 ...]
% BaselineMode is 0 since the baseline does not need to be subtracted

[FitResults,FitError] = peakfit(signal,0,0,8,1,0,10,app.first_guess,0);