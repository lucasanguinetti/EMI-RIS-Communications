clearvars -except R1 R2
close all HIDDEN
rng(10)

% This code generates Fig.6, where the SNR is plotted as a function of N
% for different non-isotropic conditions of EMI and compared with the SNR
% obtained for an algorithmically optimized RIS's phases matrix.

%Wavelength in meter
lambda = 0.1;

% Wavenumber
kappa = 2*pi/lambda;

%The width and height of an RIS element
d = lambda/4;

% Maximum mumber of elements per dimension
MaxNumOfElements = 25;

%Number of channel realizations
numOfChan = 1000;

% Vector of number of elements per dimension
NumOfElements = 1:2:MaxNumOfElements;

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
betaHd = db2pow(-inf);

% Rho value in dB - ratio between signal power and EMI power
rho = 20;

% Computing the variance of EMI from the rho values in dBm
Sigma2dBm = PowerdBm + pow2db(betaH1A/A) - rho;

% Variance of EMI in mWatt
Sigma2A = db2pow(Sigma2dBm)*A;

%% RIS Elements Coordinates
% Define the antenna elements coordinates according to a spiral
% distribution. See https://oeis.org/A174344 for reference

% Preallocate vectors
x_coordinate = zeros(1,MaxNumOfElements^2);
y_coordinate = zeros(1,MaxNumOfElements^2);

x_coordinate(1) = 0;
y_coordinate(1) = 0;

for n = 2:MaxNumOfElements^2
    x_coordinate(n) = x_coordinate(n-1) + sin(mod(floor(sqrt(4*(n-2)+1)),4)*pi/2);
    y_coordinate(n) = y_coordinate(n-1) + cos(mod(floor(sqrt(4*(n-2)+1)),4)*pi/2);
end

% Move array to the x>0,y>0 quadrant
x_coordinate = x_coordinate-min(x_coordinate);
y_coordinate = y_coordinate-min(y_coordinate);

% Scale coordinates according to the distance d
x_coordinate=x_coordinate*d;
y_coordinate=y_coordinate*d;

% Plot of the antenna coordinates
figure, plot(x_coordinate,y_coordinate,'-')
xlim([-2,MaxNumOfElements+1]*d)
ylim([-2,MaxNumOfElements+1]*d)

%% First and second channel correlation matrices
% R1: first channel

% Angle of arrival properties:
% Theta angle's mean and variance
theta_angle_mean = pi/4;
theta_angle_variance_degrees = 20;
theta_angle_variance = (theta_angle_variance_degrees*pi/180)^2;
% Phi angle's mean and variance
phi_angle_mean = 0;
phi_angle_variance = theta_angle_variance;

% Compute Correlation matrix of h1
R1 = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
    phi_angle_mean,phi_angle_variance,kappa);

% R2: second channel

% Angle of arrival properties:
% Theta angle's mean and variance
theta_angle_mean = 0;
theta_angle_variance_degrees = 20;
theta_angle_variance = (theta_angle_variance_degrees*pi/180)^2;
% Phi angle's mean and variance
phi_angle_mean = pi/4;
phi_angle_variance = theta_angle_variance;

% Compute Correlation matrix of h2
R2 = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
    phi_angle_mean,phi_angle_variance,kappa);

%% EMI correlation matrices for different variances

% Vector of EMI's variances for phi and theta
variances_vector = ([45 Inf]*pi/180).^2;

% Ready to store data
meanSNR = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_algorithm = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_rayleigh = zeros(numel(variances_vector),numel(NumOfElements));
meanSNR_noEMI = zeros(1,numel(NumOfElements));

for variances_index = 1:numel(variances_vector)
    
    variances = variances_vector(variances_index);
    
    % Rn: EMI
    
    % Angle of arrival properties:
    % Theta angle's mean and variance
    theta_angle_mean = 0;
    theta_angle_variance = variances;
    
    % Phi angle's mean and variance
    phi_angle_mean = -pi/4;
    phi_angle_variance = variances;
    
    % Generate correlation matrix of EMI
    RN = function_RCorrelated(x_coordinate,y_coordinate,theta_angle_mean,theta_angle_variance,...
        phi_angle_mean,phi_angle_variance,kappa);
    
    % Prepare to save results
    SNR = zeros(numOfChan,1);
    SNR_noEMI = zeros(numOfChan,1);
    SNR_algorithm = zeros(numOfChan,1);
    SNR_rayleigh = zeros(numOfChan,1);
    
    % Loop over number of RIS elements per dimension
    for ii = 1:numel(NumOfElements)
        
        N = NumOfElements(ii)^2;
        
        disp(['N: ',num2str(N)])
        
        % Reduced-size R1 sqrt matrix
        R1sq_i = sqrtm(R1(1:N,1:N));
        
        % Reduced-size R2 sqrt matrix
        R2sq_i = sqrtm(R2(1:N,1:N));
        
        % Reduced-size Rn matrix
        RN_i = RN(1:N,1:N);
        
        % Loop over channel realizations
        parfor kk = 1:numOfChan
            
            % Random realization of the channel defined on the paper
            h1 = sqrt(betaH1A)*R1sq_i*sqrt(.5)*(randn(N,1) + 1j*randn(N,1)); % First user channel
            h2 = sqrt(betaH2A)*R2sq_i*sqrt(.5)*(randn(N,1) + 1j*randn(N,1)); % Second user channel
            hd = sqrt(betaHd)*sqrt(.5)*(randn(1,1) + 1j*randn(1,1)); % Direct link channel            % Optimal RIS configuration against thermal noise
            
            % Optimal RIS configuration against thermal noise
            theta = diag((exp(1j*(angle(conj(h2).*h1)-angle(hd)))));
            
            % Computation of the effective channel vector g2
            g2 = theta*h2;
            
            % Current SNR value
            SNR(kk) = function_SNR(Power, SigmaW2, g2, h1, hd, Sigma2A,RN_i);
            
            % Current optimized SNR value with the iterative method and
            % upper bound
            [SNR_algorithm(kk),~,SNR_rayleigh(kk),~] = function_optimization(Power, SigmaW2, h2, h1, hd, Sigma2A, RN_i,0);
  
            if variances_index == 1
                
                SNR_noEMI(kk) = function_SNR(Power, SigmaW2, g2, h1, hd, 0, RN_i);
                
            end
        end
        
        % Save average SNR
        meanSNR(variances_index,ii) = mean(real(SNR));
        meanSNR_algorithm(variances_index,ii) = mean(real(SNR_algorithm));
        meanSNR_rayleigh(variances_index,ii) = mean(real(SNR_rayleigh));
        
        if variances_index == 1
            
            meanSNR_noEMI(1,ii) = mean(real(SNR_noEMI));
            
        end
        
    end
    
end

% Plot numerical results
% fig6 = figure(6);
% fig6.Position(3:4) = [550 350];
% hold on; box on; grid on
% set(gca,'fontsize',16);
% plot(NumOfElements.^2,smooth(meanSNR(1,:)),'LineWidth',2,'LineStyle',':','Color',[0 0 1]),
% plot(NumOfElements.^2,smooth(meanSNR_algorithm(1,:)),'LineWidth',2,'LineStyle','-.','Color',[1 0 0]),
% plot(NumOfElements.^2,smooth(meanSNR_rayleigh(1,:)),'LineWidth',2,'LineStyle','--','Color',[0 0 0]),
% plot(NumOfElements.^2,smooth(meanSNR_noEMI(1,:)),'LineWidth',2,'Color',[0 0 1])
% set(gca,'fontsize',18);
% legend({'Opt. thermal noise',...
%     'Iterative Algorithm',...
%     'Upper bound',...
%     'w/o EMI'},'Interpreter','latex','Location','best')


% Plot numerical results
fig7 = figure(7);
fig7.Position(3:4) = [550 350];
hold on; box on; grid on
set(gca,'fontsize',16);
plot(NumOfElements.^2,smooth(meanSNR(2,:)),'LineWidth',2,'LineStyle',':','Color',[0 0 1]),
plot(NumOfElements.^2,smooth(meanSNR_algorithm(2,:)),'LineWidth',2,'LineStyle','-.','Color',[1 0 0]),
plot(NumOfElements.^2,smooth(meanSNR_rayleigh(2,:)),'LineWidth',2,'LineStyle','--','Color',[0 0 0]),
plot(NumOfElements.^2,smooth(meanSNR_noEMI(1,:)),'LineWidth',2,'Color',[0 0 1])
set(gca,'fontsize',18);
legend({'Opt. thermal noise',...
    'Iterative Algorithm',...
    'Upper bound',...
    'w/o EMI'},'Interpreter','latex','Location','best')

clear Rn_i R1s R2s
