function [ R, R1_sqrt, R2_sqrt ] = function_CorrMatComputation_Iso(sqrtN, d, lambda, betaH1A, betaH2A)

%% RIS Elements Coordinates
% Define the elements coordinates according to a square lattice
gridPoints = (0:sqrtN-1)*d;

[X,Y] = meshgrid(gridPoints,gridPoints); % coordinates expressed as matrices

locations = X(:)+1i*Y(:); % coordinates expressed as complex vector


%% Spatial Correlation Matrices

% Isotropically irradiated RIS correlation matrix
R = sinc(2*abs(locations - transpose(locations))/lambda);
R_sqrt = sqrtm(R);

% User 1 channel correlation matrix times betaH1A
% R1 = betaH1A*R;
R1_sqrt = sqrt(betaH1A)*R_sqrt; %sqrtm(R1);

% User 2 channel correlation matrix times betaH2A
% R2 = betaH2A*R;
R2_sqrt = sqrt(betaH2A)*R_sqrt; %sqrtm(R2);
