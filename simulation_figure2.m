close all;
clear;

% This code generates Fig.2, where the SNR is plotted as a function of N
% for different values of rho

%Wavelength in meter
lambda = 0.1;

% Wavenumber
kappa = 2*pi/lambda;

%The width and height of an RIS element
d = lambda/4;

% Maximum mumber of elements per dimension
MaxNumOfElements = 50;

%Number of channel realizations
numOfChan = 1000;

% Vector of number of elements per dimension
NumOfElements = 10:5:MaxNumOfElements;

%% System Parameters - as defined in the letter

% Bandidth
Bandwidth = 1e6;

% RIS element area
A = d.^2;

% Total Radiated Power in dBm
PowerdBm = 23;

% Total Radiated Power in mWatt
Power = db2pow(PowerdBm);

% Thermal noise in dBm
SigmaW2dBm = pow2db(Bandwidth)-174; % -114 dBm;

% Thermal noise in Watt (sigma2_w)
SigmaW2 = db2pow(SigmaW2dBm);

% SNR
Psigma2 = Power/SigmaW2;

% SNR in dB
Psigma2dB = pow2db(Psigma2);

% Channel gain h1
betaH1A = db2pow(-48)*A;

% Channel gain h2
betaH2A = db2pow(-38)*A;

% Channel gain direct link
betaHd = db2pow(-inf);

% Vector containing rho values - ratio between signal power and EMI power
vect_rho = [10 20 30];

% Ready to store data for w/o and w/ EMI cases
meanSNR = zeros(numel(vect_rho), numel(NumOfElements));
meanSNR_noEMI = zeros(numel(vect_rho), numel(NumOfElements));

% Loop over RIS elements
for ii  = 1:numel(NumOfElements)
    
    sqrtN = NumOfElements(ii);
    
    N = sqrtN^2;
    
    disp(['Sqrt N: ',num2str(sqrtN)])
    
    % Generate correlation matrices
    [ Rn, R1_sqrt, R2_sqrt ] = function_CorrMatComputation_Iso(sqrtN, d, lambda, betaH1A, betaH2A);
     
    % Loop over rho values
    for index = 1:numel(vect_rho)
        
        % Current value of rho
        rho = vect_rho(index);
        
        % Computing the variance of EMI from the rho values in dBm
        Sigma2dBm = PowerdBm + pow2db(betaH1A/A) - rho;
        
        % Variance of EMI in mWatt
        Sigma2A = db2pow(Sigma2dBm)*A;

        % Computation of SNR for different channel realizations
        [SNR_noEMI, SNR,~,~] = function_SNR_v2(N, numOfChan,...
            Power, SigmaW2, Sigma2A, Rn, R1_sqrt, R2_sqrt, betaHd);
        
        % Average over channel realizations
        meanSNR(index,ii) = mean(real(SNR));
        meanSNR_noEMI(index,ii) = mean(real(SNR_noEMI));
        
    end
    
end

% Plot numerical results
fig2 = figure(2);hold on; box on; grid on;
plot(NumOfElements,(meanSNR_noEMI(1,:)),'k-','LineWidth',3);
plot(NumOfElements,(meanSNR(3,:)),'r-.','LineWidth',3);
plot(NumOfElements,(meanSNR(2,:)),'r-.','LineWidth',3);
plot(NumOfElements,(meanSNR(1,:)),'r-.','LineWidth',3);
set(gca,'fontsize',18);
set(gca,'Yscale','log');
xlabel('Number of elements (per dimension)','Interpreter','latex');
ylabel('Average SNR','Interpreter','latex');
legend({'w/o EMI','w/ EMI'},'Interpreter','latex','Location','NorthWest');
xlim([10 MaxNumOfElements])
legend
fig2.Position(3:4) = [550 350];

clear Rn R1_sqrt R2_sqrt

%save matlab_fig2.mat
