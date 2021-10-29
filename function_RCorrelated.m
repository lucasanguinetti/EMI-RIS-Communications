function [RMatrix] = function_RCorrelated(X,Y,t_mean,t_var,p_mean,p_var,k)



% PDF function as in Eq.XX
pdf_function = @(theta,phi)  exp(-(pi-abs(pi-abs(t_mean-theta))).^2/t_var).*exp(-(p_mean-phi).^2/p_var).*(cos(theta)./(2.*pi));

% ACF function for user EMI
acf_function = @(y,z) integral2(@(theta,phi) pdf_function(theta,phi).*exp(-1j.*k.*(sin(phi).*cos(theta).*y + sin(theta).*z)),...
    -pi/2,pi/2,-pi/2,pi/2);

acf_function_normalized = @(y,z) acf_function(y,z)./acf_function(0,0);

% Correlation matrix
% The antennas coordinates need to be expressed as vectors
RMatrix = RMatrixWaitbar(X(:),Y(:),acf_function_normalized,'R matrix');



end

%% Function R
function R = RMatrixWaitbar(X,Y,acf_fun,stage)

% This function compute the R matrix entries corrensonding to the antenna
% points defined by the coordinates vectors X and Y under an
% autocorrelation functione defined in acf_fun. stage is used to apply a
% label on the waitbar

D = parallel.pool.DataQueue;
h = waitbar(0, stage);
afterEach(D, @nUpdateWaitbar);

p = 1;

% Coordinates under matricial form
X_mat = X-X';
Y_mat = Y-Y';

% Empty matrix R ready to store data
R = X_mat*0;
N = size(X_mat,1);
parfor m = 1:N
    for n = 1:N
        R(m,n) = acf_fun(X_mat(m,n),Y_mat(m,n));
        send(D,m)
    end
    
end
clear X_mat Y_mat

    function nUpdateWaitbar(~)
        waitbar(p/(N*N), h);
        p = p + 1;
    end
close(h)
end


