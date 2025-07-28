close all;
clear;
clc;

% Select the range of SNR values
SNRdB = -10;
SNR = db2pow(SNRdB);

% Select range of the number of antennas
Mall = 8:1:128;
K = 8;
NTry = 1000;

% Prepare to save simulation results
sumrate_LMMSE = zeros(length(Mall),1);
sumrate_MRC = zeros(length(Mall),1);
sumrate_ZF = zeros(length(Mall),1);

% Temporary variables for Monte Carlo iterations
sumrate_LMMSE_temp = zeros(length(Mall),NTry);
sumrate_MRC_temp = zeros(length(Mall),NTry);
sumrate_ZF_temp = zeros(length(Mall),NTry);

% Store all Monte Carlo results for individual users
rate_LMMSE = zeros(K,length(Mall),NTry);
SINR_LMMSE = zeros(K,length(Mall),NTry);
rate_MRC = zeros(K,length(Mall),NTry);
SINR_MRC = zeros(K,length(Mall),NTry);
SINR_ZF = zeros(K,length(Mall),NTry);
rate_ZF = zeros(K,length(Mall),NTry);

for m = 1:length(Mall)
    M = Mall(m);
    
    % Monte Carlo simulation
    for n = 1:NTry
        % Generate random channel matrix for each iteration
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);

        powerCoefNorm = ones(K,1);
        
        % LMMSE receiver
        for k = 1:K
            SINR_LMMSE(k,m,n) = SNR*powerCoefNorm(k)*real(H(:,k)'*((SNR*(H*diag(powerCoefNorm)*H'-powerCoefNorm(k)*H(:,k)*H(:,k)')+eye(M))\H(:,k)));
            rate_LMMSE(k,m,n) = log2(1+ SINR_LMMSE(k,m,n));
        end
        sumrate_LMMSE_temp(m,n) = sum(rate_LMMSE(:,m,n));
        
        % MRC receiver
        W_MRC = H;
        for k = 1:K
            SINR_MRC(k,m,n) = SNR*powerCoefNorm(k)*abs(H(:,k)'*W_MRC(:,k))^2/real(W_MRC(:,k)'*(SNR*(H*diag(powerCoefNorm)*H'-powerCoefNorm(k)*H(:,k)*H(:,k)')+eye(M))*W_MRC(:,k));
            rate_MRC(k,m,n) = log2(1+ SINR_MRC(k,m,n));
        end
        sumrate_MRC_temp(m,n) = sum(rate_MRC(:,m,n));
        
        % ZF receiver
        
        W_ZF = H/(H'*H);
        for k = 1:K
            SINR_ZF(k,m,n) = SNR * powerCoefNorm(k) * abs(H(:,k)' * W_ZF(:,k))^2 / real(W_ZF(:,k)' * (SNR * (H * diag(powerCoefNorm) * H' - powerCoefNorm(k) * H(:,k) * H(:,k)') + eye(M)) * W_ZF(:,k));
            rate_ZF(k,m,n) = log2(1+ SINR_ZF(k,m,n));
        end
        sumrate_ZF_temp(m,n) = sum(rate_ZF(:,m,n));
            
    end
    
    % Calculate average over Monte Carlo iterations
    sumrate_LMMSE(m) = mean(sumrate_LMMSE_temp(m,:));
    sumrate_MRC(m) = mean(sumrate_MRC_temp(m,:));
    sumrate_ZF(m) = mean(sumrate_ZF_temp(m,:));
    fprintf('Completed M = %d antennas (%d/%d)\n', M, m, length(Mall));   
end

% Plot simulation results
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure;
hold on; box on; grid on;
plot(Mall, sumrate_LMMSE, 'k--', 'LineWidth', 2);
plot(Mall, sumrate_MRC, 'r-o', 'LineWidth', 2);
plot(Mall, sumrate_ZF, 'b-', 'LineWidth', 2);
xlabel('Number of antennas $(M)$','Interpreter','latex');
ylabel('Sum rate [bit/symbol]','Interpreter','latex');
set(gca,'fontsize',16);
legend({'LMMSE','MRC','ZF'},'Interpreter','latex','Location','NorthWest');
title('Hiệu suất các kỹ thuật kết hợp thu trong Massive MIMO đa người dùng đường lên','Interpreter','latex');

grid on;
xlim([Mall(1) Mall(end)])