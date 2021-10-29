function [SNR] = function_SNR(Power, SigmaW2, g2, h1, hd, Sigma2A,Rn)

% This function computes the istantaneous SNR receiving a input the following
% Power: the scalar value of the transmitted power
% SigmaW2: the scalar value of the thermal noise power
% g2: the second user channel, usually accounting for the RIS phase shift matrix
% h1: the first user channel
% Sigma2A: the scalar value of the EMI magnitude times the RIS element area
% Rn: the correlation matrix of the EMI


SNR = Power*abs(g2'*h1 + hd)^2/(Sigma2A*(g2'*Rn*g2)+SigmaW2);

end

