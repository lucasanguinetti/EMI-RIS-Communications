function [SNR_nonoise_vett, SNR_noise_vett,SNR_nonoise_vett_THETA, SNR_noise_vett_THETA] = function_SNR_v2(N,K,Power,SigmaW2,Sigma2A,R,Rh1_sqrt,Rh2_sqrt,betaHd)


% This function computes K times the istantaneous SNR receiving a input the following
% N: the number of array elements
% K: the number of channel realizations
% Power: the scalar value of the transmitted power
% SigmaW2: the scalar value of the thermal noise power
% Sigma2A: the scalar value of the EMI magnitude times the RIS element area
% Rn: the correlation matrix of the EMI
% Rh1_sqrt: the first user channel root correlation matrix
% Rh2_sqrt: the second user channel root correlation matrix
% betaHd: the average magnitude of the direct link amplitude

SNR_nonoise_vett = zeros(K,1);
SNR_noise_vett = zeros(K,1);
SNR_nonoise_vett_THETA = zeros(K,1);
SNR_noise_vett_THETA = zeros(K,1);

for k = 1:K
    h1 = sqrt(.5)*(Rh1_sqrt)*(randn(N,1) + 1j*randn(N,1));
    h2 = sqrt(.5)*(Rh2_sqrt)*(randn(N,1) + 1j*randn(N,1));
    hd = sqrt(.5)*(sqrt(betaHd))*(randn(1,1) + 1j*randn(1,1));
    
    
  
    %% SNRs
    % Optimized phases SNR computation
    SNR_nonoise_vett(k) = Power*(sum(abs(h1.*h2),1)+abs(hd)).^2./(SigmaW2);
    SNR_noise_vett(k) = Power*(sum(abs(h1.*h2),1)+abs(hd)).^2./(Sigma2A*(h2'*R*h2) + SigmaW2); %Power*sum(abs(h.*g)+abs(hd),1).^2./((g'*betaN*g) + Sigma2);
    
    % Explicit phases optimization through theta (results are as above)
    theta = exp(1i*(angle(h2.*conj(h1)) - angle(hd)));
    g2 = diag(theta')*(h2);
    
    SNR_nonoise_vett_THETA(k) = Power*abs(g2'*h1 + hd).^2./(SigmaW2);
    SNR_noise_vett_THETA(k) = Power*abs(g2'*h1 + hd).^2./(Sigma2A*(g2'*R*g2) + SigmaW2); %Power*sum(abs(h.*g)+abs(hd),1).^2./((g'*betaN*g) + Sigma2);
    
    
end

end

