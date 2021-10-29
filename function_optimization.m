function [SNR_optimized_output,SNR_optimized_thermal_noise_output,SNR_rayleigh_quotient_maximizer_output,SNR_rayleigh_quotient_maximizer_phase_output] = ...
    function_optimization(Power, SigmaW2, h2, h1, hd, Sigma2A,Rn,show)

%% Optimization algorithm loop

Number_total_elements = numel(h1);
Number_elements_per_side = sqrt(Number_total_elements);

% Auxiliary matrices as defined in the paper
H2 = diag(h2);
A = H2'*(h1*h1')*H2;
B = H2'*Rn*H2;
C = Sigma2A/SigmaW2*B + (1/Number_total_elements)*eye(Number_total_elements);
D = sqrtm(inv(C'))*A*sqrtm(inv(C));
Dsq = sqrtm(D); % D^(1/2)
alpha = .5/norm(Dsq'*Dsq,2);

% Initialization
theta_optimized = ones(Number_total_elements,1);
iterations_number = 50;

% Ready to store data
SNR_optimized = zeros(1,iterations_number);
SNR_optimized_thermal_noise = zeros(1,iterations_number);
SNR_rayleigh_quotient_maximizer = zeros(1,iterations_number);
SNR_rayleigh_quotient_maximizer_phase = zeros(1,iterations_number);


% Loop
for i = 1:iterations_number
    
%     disp(['Iteration: ',num2str(i)])
    
    % Meaning of Subscripts:
    % _optimized: values resulting after Gradient Projrection optimization.
    % _optimized_thermal_noise: values resulting after the optimization of
    % theta in the absence of EMI.
    % _rayleigh_quotient_maximizer: values resulting after the maximization
    % of the Rayleigh quatient in the notes.
    % _rayleigh_quotient_maximizer_phase: values computed with the previous
    % case theta's phase.
    
    % zeta and theta_optimized are the ones defined in the note
    zeta = sqrtm((C))*theta_optimized + alpha*(Dsq'*Dsq)*sqrtm((C))*theta_optimized;
    theta_optimized = exp(1j*angle(sqrtm(inv(C))*zeta));
    
    % Theta vector, as defined in the paper
    theta_optimized_thermal_noise = exp(1j*(angle(conj(h2).*(h1)))); % original
    theta_rayleigh_quotient_maximizer = ((C\(H2'*h1)));
    theta_rayleigh_quotient_maximizer = Number_elements_per_side*theta_rayleigh_quotient_maximizer/norm(theta_rayleigh_quotient_maximizer);
    theta_rayleigh_quotient_maximizer_phase = exp(1j*angle(theta_rayleigh_quotient_maximizer));
    
    % g2 = diag(theta)*h2, as defined in the paper
    g2_optimized = diag(theta_optimized)*h2;
    g2_optimized_thermal_noise = diag(theta_optimized_thermal_noise)*h2;
    g2_rayleigh_quotient_maximizer = diag(theta_rayleigh_quotient_maximizer)*h2;
    g2_rayleigh_quotient_maximizer_phase = diag(theta_rayleigh_quotient_maximizer_phase)*h2;
    
    
    % SNR = (Power/SigmaW2)*abs(g2'*h1 + hd)^2/(gamma*(g2'*Rn*g2)+1); with
    % gamma = Sigma2A/SigmaW2;
    SNR_optimized(i) = function_SNR(Power, SigmaW2, g2_optimized, h1, hd, Sigma2A,Rn);
    SNR_optimized_thermal_noise(i) = function_SNR(Power, SigmaW2, g2_optimized_thermal_noise, h1, hd, Sigma2A,Rn);
    SNR_rayleigh_quotient_maximizer(i) = function_SNR(Power, SigmaW2, g2_rayleigh_quotient_maximizer, h1, hd, Sigma2A,Rn);
    SNR_rayleigh_quotient_maximizer_phase(i) = function_SNR(Power, SigmaW2, g2_rayleigh_quotient_maximizer_phase, h1, hd, Sigma2A,Rn);
    
    
    
end

% The last element of the optimization iterations are delivered as output
% of the function
SNR_optimized_output = SNR_optimized(end);
SNR_optimized_thermal_noise_output = SNR_optimized_thermal_noise(end);
SNR_rayleigh_quotient_maximizer_output = SNR_rayleigh_quotient_maximizer(end);
SNR_rayleigh_quotient_maximizer_phase_output = SNR_rayleigh_quotient_maximizer_phase(end);


if show
    f = figure;
    hold on, grid on
    plot(1:numel(SNR_optimized),SNR_optimized,'LineWidth',2)
    plot(1:numel(SNR_optimized_thermal_noise),SNR_optimized_thermal_noise,'LineWidth',2)
    plot(1:numel(SNR_rayleigh_quotient_maximizer),SNR_rayleigh_quotient_maximizer,'LineWidth',2)
    plot(1:numel(SNR_rayleigh_quotient_maximizer_phase),SNR_rayleigh_quotient_maximizer_phase,'LineWidth',2)
    legend({'Optimized','Optimized for thermal noise','Rayleigh quotient maximizer','Quotient maximizer phase'},'Location','best','Interpreter','latex')
    ylabel('SNR','Interpreter','latex')
    f.Position(3:4) = [550, 368];
end

end

