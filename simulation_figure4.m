close all;
clear;

% This code generates Fig. 4, in which we compute the asymptotic value of
% SNR

% Wavelength in meter
lambda = 0.1;

% Wavenumber
kappa = 2*pi/lambda;

%The width and height of an RIS element
d = lambda/4;

% Maximum mumber of elements per dimension
MaxNumOfElements = 100;

%Number of channel realizations
numOfChan = 10;

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

% Channel gain of direct link
vector_betaHd = [-inf -100  -80];

% Rho value in dB - ratio between signal power and EMI power
rho = 20;

% Computing the variance of EMI from the rho values in dBm
Sigma2dBm = PowerdBm + pow2db(betaH1A/A) - rho;

% Variance of EMI in mWatt
Sigma2A = db2pow(Sigma2dBm)*A;

% Ready to store data
SNR = zeros(numel(vector_betaHd), numel(NumOfElements));

% Loop over RIS elements
for ii  = 1:numel(NumOfElements)
    
    sqrtN = NumOfElements(ii);
    
    N = sqrtN^2;
    
    disp(['Sqrt N: ',num2str(sqrtN)])
    
    % Generate correlation matrices
    [ Rn, R1_sqrt, R2_sqrt ] = function_CorrMatComputation_Iso(sqrtN, d, lambda, betaH1A, betaH2A);
    
    % Loop over direct link channel gains
    for index = 1:numel(vector_betaHd)
        
        betaHd = db2pow(vector_betaHd(index)); % Direct link channel

        % Ready to store instantenous SNR values
        SNR_vect = zeros(numOfChan,1);
        
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
            SNR_vect(kk) = function_SNR(Power, SigmaW2, g2, h1, hd, Sigma2A,Rn);
            
        end
        
        % Compute the average SNR
        SNR(index, ii) = mean(real(SNR_vect));
        
    end
    
end


% %% Computation of Asymptotic Value
%
% % Asymptotic line
% NAsympSq = 200; % Large values on N. 100 is given by the MacBook RAM limit
% % Redefinition of the antennas lattice and correlation matrix as in
% % previous line. Here, the purpose is to estimate the asympotic
% % limit, for N->inf, of Eq.XX
% gridPoints = (0:NAsympSq-1)*d;
% [X,Y] = meshgrid(gridPoints,gridPoints);
% locations = X(:)+1i*Y(:);
% NAsymp = length(locations);
% sig = Sigma2A/A;
% betaH1 = betaH1A/A;
% R = sinc(2*abs(locations - transpose(locations))/lambda);
% Alpha = (1/NAsymp)* trace(R*R);
% Asymp = ones(1,numel(NumOfElements))*Power*betaH1*(pi/4).^2/sig/Alpha; % Froma Eq.XX on the paper
% clear R


Asymp = ones(1,numel(NumOfElements))*(10*log10(SNR(1,end)./(NumOfElements(end).^2)));

fig4 = figure(4);
hold on; box on; grid on;
plot(NumOfElements',(10*log10(SNR(1,:)./(NumOfElements.^2))),'r-.','LineWidth',3); % N_ant is the square root of total N
plot(NumOfElements',(10*log10(SNR(2,:)./(NumOfElements.^2))),'b--','LineWidth',3);
plot(NumOfElements',(10*log10(SNR(3,:)./(NumOfElements.^2))),'k:','LineWidth',3);
plot(NumOfElements',Asymp,'k-','LineWidth',3);
set(gca,'fontsize',18);
set(gca,'Yscale','lin');
xlabel('Number of elements (per dimension)','Interpreter','latex');
ylabel('Normalized SNR [dB]','Interpreter','latex');
legend({'EM-Noise ','EM-Noise bound'},'Interpreter','latex','Location','Best');
xlim([10 100]);
ylim([-0 40]);
fig4.Position(3:4) = [550 350];
legend({'$\beta_{\rm d} = -\infty$',  '$\beta_{\rm d} = -100$ dB', '$\beta_{\rm d} = -80$ dB'},'Interpreter','latex','FontSize',18);

clear Rn R1_sqrt R2_sqrt

%save matlab_fig4.mat
