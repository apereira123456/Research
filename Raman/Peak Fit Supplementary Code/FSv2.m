clear
clc
close all
noise=0.001;
heatingRate = 10;
%% READ DATA
% Specify the folder where the files live.
%myFolder = 'C:\Users\Eoin\Dropbox\Experiments\massspec\Composition Variation\Plot';
myFolder = pwd;
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf(...
        'Error: The following folder does not exist:\n%s\nPlease specify a new folder.',...
        myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
        % User clicked Cancel
        return;
    end
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.xlsx'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
number_of_files = length(theFiles);
ramp = cell(number_of_files, 1);
filenames = cell(number_of_files, 1);
for i = 1 : length(theFiles)
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(theFiles(i).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    filenames{i} = baseFileName(1:end-5);
    filenames{i} = strrep(filenames{i},'_',' ');
    numVars = 4;
    varTypes = {'double','double','double','double'} ;
    varNames = {'Time','Temperature','HeatFlow','Weight'} ;
    opts = spreadsheetImportOptions('NumVariables',numVars,...
        'VariableNames',varNames,...
        'VariableTypes',varTypes,...
        'DataRange','A4');
    ramp{i} = readtable(fullFileName,opts,'Sheet', 3);
    
end


isocat = cell(number_of_files,1);
lag = 1000;
isocat_short = cell(number_of_files,1);
for i = 1 : length(theFiles)
    ramp{i}.NormWeight = (ramp{i}.Weight - min(ramp{i}.Weight))./(max(ramp{i}.Weight)-min(ramp{i}.Weight));
    isocat{i} = ramp{i};
    isocatstore = isocat{i};
    isocat_short{i} = isocatstore(lag:lag:end, :);
    isocat_short{i}.Deriv = gradient(isocat_short{i}.Weight)./gradient(isocat_short{i}.Temperature);
    deriv = isocat_short{i}.Deriv;
    %isocat_short{i}.Deriv(1:150) = sgolayfilt(isocat_short{i}.Deriv(1:150),2,11);
    %isocat_short{i}.Deriv = sgolayfilt(isocat_short{i}.Deriv,2,21);
    %ramp{i}.NormWeight = (ramp{i}.Weight - min(ramp{i}.Weight))./(max(ramp{i}.Weight)-min(ramp{i}.Weight));
end
%% Assign data to Variables
M100temp = isocat_short{1}.Temperature+273.15;
M100deriv = isocat_short{1}.Deriv;
global C0
M100weight = isocat_short{1}.NormWeight;
M100weight(1);
C0 = M100weight(1);
%% PeakFitting
global PEAKHEIGHTS AUTOZERO
format short g
format compact
warning off all
NumPeaks=3;
AUTOZERO=0;
n=length(M100temp );
% Starting first guess with small random variaiton for each repeat
%newstart=[607 55 0 670 45 -.15 755 60 -.3 500 30 -.25];
newstart=[607 55 0 670 45 -.15 755 60 -.3];
for parameter=1:3:3*NumPeaks
    newstart(parameter)=newstart(parameter)*(1+randn/100);
    newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
end

options = optimset('TolX',.001,'Display','off','MaxFunEvals',2000);
tic
TrialParameters=fminsearch(@(lambda)(fitFS(lambda,M100temp,M100deriv)),newstart);
toc
disp('       Center     Width       shape      center     width     shape')
disp(TrialParameters)
alphas=[TrialParameters(3) TrialParameters(6)]
% Construct model from Trial parameters
A=zeros(NumPeaks,n);
dcdtfitclean=zeros(NumPeaks,n);
dcdtfit=zeros(NumPeaks,n);

for m=1:NumPeaks
    A(m,:)=FS(M100temp,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
    dcdtfitclean(m,:) = PEAKHEIGHTS(m).*A(m,:);
    dcdtfit(m,:) = dcdtfitclean(m,:)+(noise.*randn(size(M100temp)))';
end
% Multiplies each row by the corresponding amplitude and adds them up
if AUTOZERO==3
    baseline=PEAKHEIGHTS(1);
    Heights=PEAKHEIGHTS(2:1+NumPeaks);
    model=Heights'*A+baseline;
else
    model=PEAKHEIGHTS'*A;
    Heights=PEAKHEIGHTS;
    baseline=0;
end
% Compare trial model to data segment and compute the fit error
MeanFitError=100*norm(M100deriv-model)./(sqrt(n)*max(M100deriv))
%plot(M100temp,M100deriv,'.k',M100temp,model,'r')
%% Plot Peak Fit
% figure()
% hold on
% plot(M100temp,M100deriv,'LineWidth',2)
% for i = 1:NumPeaks
%     plot(M100temp,PEAKHEIGHTS(i).*A(i,:),'LineStyle','--','LineWidth',1.25);
% end
% ax = gca;
% ax.FontSize = 15;
% xlabel('Temperature (K)')
% ylabel('Rate of Mass Loss (%/K)')
%%
close all
heatingRate = 1;
Beta = heatingRate/60;
global INVBeta
INVBeta = 1/Beta;
Tspan = M100temp(50:400)'+273.15;
dcdtraw = dcdtfit(1,50:400);
% param0 = [5 5000 1];
% param = cell(1,1);
mag = [1 10 100 1000 10000 100000]
mag2 = [1 5 10 1000 10000 100000]
%rand*mag2(randsample(3,1))
params = [rand*mag(randsample(4,1)) rand*mag(randsample(4,1))];
param0 = params;
%param0 = [randsample(125,1) randsample(30000,1)]
resid = cell(1,1);
J = cell(1,1);
CovB = cell(1,1);
opts = statset('MaxIter',10000, 'Display', 'iter', 'TolX', eps, 'TolFun', 1e-8);% Options: displays iteration #, Max iteration 1000, the tolerance in the x direction is nul.
%myfunc(param0,Tspan)
[param{1}, resid{1}, J{1}, CovB{1}] = nlinfit(Tspan,dcdtraw,@myfunc,param0,opts);
figure()
hold on
plot(Tspan,  myfunc(param{1},Tspan))

%plot(Tspan, dcdtraw)


%% Functions
function err = fitFS(lambda, temp, deriv)
global PEAKHEIGHTS AUTOZERO
A = zeros(length(temp), round(length(lambda)/3));
for j = 1:length(lambda)/3
    A(:,j) = FS(temp,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
end
if AUTOZERO==3,A=[ones(size(deriv))' A];end
PEAKHEIGHTS = A\deriv;
z = A*PEAKHEIGHTS;
err = norm(z-deriv);
end
function g=FS(x,center,width,shape)
g = exp(-log(2).*(log(1+2.*shape.*((x-center)./width))./shape).^2);
g(g~=real(g))=0;
end
function output = myfunc(param,Tspan)
global C0
global INVBeta
Ca_new = zeros(length(Tspan),1);
display(sprintf('%1.3f, %2.0f, %3.3f', param))% %4.3f',param))
% we just rename the variable here for use in the ODE
lnA = param(1);
ER = param(2);
%n = param(3);
n =1/3;
%Ca0 = param(4);
Ca_new(1) = C0;
%Ca = Ca0;


dcdt = zeros(length(Tspan),1);
%[T,Ca] = ode45(@ode,Tspan,Ca0);
for i = 1:(length(Tspan)-1)
    if Ca_new(i) < 0
        Ca_new(i) = 0;
    end %if
    %dcdt(i) = -INVBeta.*exp(lnA-(ER./Tspan(i))).*Ca_new(i).^(n);
    %dcdt(i) = -INVBeta.*exp(lnA-(ER./Tspan(i))).*(Ca_new(i).*((-log(Ca_new(i))).^(1-(1/n))));
    %dcdt(i) = -INVBeta.*exp(lnA-(ER./Tspan(i))).*1./(-log(Ca_new(i)));
    %dcdt(i) = -INVBeta.*exp(lnA-(ER./Tspan(i))).*((3/2).*(Ca_new(i).^(n)-1));
    Ca_new(i)
    Ca_new(i+1) = Ca_new(i) + dcdt(i) .* (Tspan(i+1)-Tspan(i));
end %for
output = dcdt';
end %myfunc