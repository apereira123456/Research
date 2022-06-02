function baseline = baselinefit(signal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:    signal (baseline plus peaks)
% Output:   baseline
% This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% Copyright 2013 Semih Agcaer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poly_degree = 6;
fit = polyfit((1:numel(signal)),signal, poly_degree);          
baseline_temp =(polyval(fit,1:numel(signal)));
signal_new = signal;
signal_new(signal_new>baseline_temp) = baseline_temp((signal_new>baseline_temp));
poly_degree = 11;
fit = polyfit((1:numel(signal)),signal_new, poly_degree);          
baseline =(polyval(fit,1:numel(signal_new)));
