close all;
clear;

% This code generates Fig.3, where the SNR is plotted as a function of N
% for different values of direct link channel gain

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
NumOfElements = 10:2:MaxNumOfElements;

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

% Channel gain of direct link
vector_betaHd = [-inf -100  -80];

% Rho value in dB - ratio between signal power and EMI power
rho = 20;

% Computing the variance of EMI from the rho values in dBm
Sigma2dBm = PowerdBm + pow2db(betaH1A/A) - rho;

% Variance of EMI in mWatt
Sigma2A = db2pow(Sigma2dBm)*A;


% Ready to store data for w/o RIS, w/o EMI and w/ EMI cases
meanSNR = zeros(numel(vector_betaHd), numel(NumOfElements));
meanSNR_noRIS = zeros(numel(vector_betaHd), numel(NumOfElements));
meanSNR_noEMI = zeros(numel(vector_betaHd), numel(NumOfElements));

% Loop over RIS elements
for ii  = 1:numel(NumOfElements)
    
    sqrtN = NumOfElements(ii);
    
    N = sqrtN^2;
    
    disp(['Sqrt N: ',num2str(sqrtN)])
    
    % Generate correlation matrices
    [ Rn, R1_sqrt, R2_sqrt ] = function_CorrMatComputation_Iso(sqrtN, d, lambda, betaH1A, betaH2A);
    
    % Loop over direct link channel gains
    for index = 1:numel(vector_betaHd)
        
        % Current channel gain of direct link in dB
        betaHddB = vector_betaHd(index);
        
        % Current channel gain of direct link in mWatt
        betaHd = db2pow(betaHddB);

        % Ready to store instantenous SNR values
        SNR = zeros(numOfChan,1);
        SNR_noRIS = zeros(numOfChan,1);
        SNR_noEMI = zeros(numOfChan,1);
       
        % Loop over channel realizations
        for kk = 1:numOfChan
            
            % Generate Channel Vectors
            h1 = R1_sqrt*sqrt(.5)*(randn(N,1) + 1j*randn(N,1));
            h2 = R2_sqrt*sqrt(.5)*(randn(N,1) + 1j*randn(N,1));
            hd = sqrt(betaHd)*sqrt(.5)*(randn(1,1) + 1j*randn(1,1));
            
            % Optimal RIS configuration against thermal noise
            theta = diag((exp(1j*(angle(conj(h2).*h1)-angle(hd)))));
            
            % Computation of the effective channel vector g2
            g2 = theta*h2;
            
            % Compute the SNR value for the current channel realization
            SNR(kk) = function_SNR(Power, SigmaW2, g2, h1, hd, Sigma2A,Rn);
            
            % Compute the SNR value w/o RIS for the current channel realization
            SNR_noRIS(kk) = function_SNR(Power, SigmaW2, 0*g2, 0*h1, hd, 0,Rn);
            
            % Compute the SNR value w/o EMI for the current channel realization
            SNR_noEMI(kk) = function_SNR(Power, SigmaW2, g2, h1, hd, 0,Rn);
            
        end
        
        % Average over channel realizations
        meanSNR(index, ii) = mean(real(SNR));
        meanSNR_noRIS(index, ii) = mean(real(SNR_noRIS));
        meanSNR_noEMI(index, ii) = mean(real(SNR_noEMI));
        
    end
    
end

% Plot numerical results
fig3 = figure(3);
for index = 2:numel(vector_betaHd)
    hold on; box on; grid on;
    plot(NumOfElements,smooth(meanSNR_noEMI(index,:),7),'k-','LineWidth',3);
    plot(NumOfElements,smooth(meanSNR(index,:),7),'r-.','LineWidth',3);
    plot(NumOfElements,ones(1,numel(NumOfElements))*(meanSNR_noRIS(index,1)),'b:','LineWidth',3);
    set(gca,'fontsize',18);
    set(gca,'Yscale','log');
    xlabel('Number of elements (per dimension)','Interpreter','latex');
    ylabel('Average SNR','Interpreter','latex');
end
legend({'w/o EMI','w/ EMI','w/o RIS'},'Interpreter','latex','Location','Best');
xlim([10 MaxNumOfElements])
ylim([1e3 1e7])
fig3.Position(3:4) = [550 350];

clear Rn R1_sqrt R2_sqrt

%save matlab_fig3.mat
