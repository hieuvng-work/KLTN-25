%% khao sat SNRdB với mô hình kênh Rayleigh fading
close all;
clear;
clc;

% Select the range of SNR values
SNRdB = -10:1:30;
SNR = db2pow(SNRdB);

% Select range of the number of antennas
M = 64;
K = 8;
NTry = 1e4;

% Prepare to save simulation results
sumrate_LMMSE = zeros(length(SNR),1);
sumrate_MRC = zeros(length(SNR),1);
sumrate_ZF = zeros(length(SNR),1);

% Temporary variables for Monte Carlo iterations
sumrate_LMMSE_temp = zeros(length(SNR),NTry);
sumrate_MRC_temp = zeros(length(SNR),NTry);
sumrate_ZF_temp = zeros(length(SNR),NTry);

% Store all Monte Carlo results for individual users
rate_LMMSE = zeros(K,length(SNR),NTry);
SINR_LMMSE = zeros(K,length(SNR),NTry);
rate_MRC = zeros(K,length(SNR),NTry);
SINR_MRC = zeros(K,length(SNR),NTry);
SINR_ZF = zeros(K,length(SNR),NTry);
rate_ZF = zeros(K,length(SNR),NTry);

for s = 1:length(SNRdB)
    % Monte Carlo simulation
    for i = 1:NTry
        % Generate random channel matrix for each iteration
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);

        powerCoefNorm = ones(K,1);

        % LMMSE receiver
        for k = 1:K
            SINR_LMMSE(k,s,i) = SNR(s)*powerCoefNorm(k)*real(H(:,k)'*((SNR(s)*(H*diag(powerCoefNorm)*H'-powerCoefNorm(k)*H(:,k)*H(:,k)')+eye(M))\H(:,k)));
            rate_LMMSE(k,s,i) = log2(1+ SINR_LMMSE(k,s,i));
        end
        sumrate_LMMSE_temp(s,i) = sum(rate_LMMSE(:,s,i));

        % MRC receiver
        W_MRC = H;
        for k = 1:K
            SINR_MRC(k,s,i) = SNR(s)*powerCoefNorm(k)*abs(H(:,k)'*W_MRC(:,k))^2/real(W_MRC(:,k)'*(SNR(s)*(H*diag(powerCoefNorm)*H'-powerCoefNorm(k)*H(:,k)*H(:,k)')+eye(M))*W_MRC(:,k));
            rate_MRC(k,s,i) = log2(1+ SINR_MRC(k,s,i));
        end
        sumrate_MRC_temp(s,i) = sum(rate_MRC(:,s,i));

        % ZF receiver

        W_ZF = H/(H'*H);
        for k = 1:K
            SINR_ZF(k,s,i) = SNR(s) * powerCoefNorm(k) * abs(H(:,k)' * W_ZF(:,k))^2 / real(W_ZF(:,k)' * (SNR(s) * (H * diag(powerCoefNorm) * H' - powerCoefNorm(k) * H(:,k) * H(:,k)') + eye(M)) * W_ZF(:,k));
            rate_ZF(k,s,i) = log2(1+ SINR_ZF(k,s,i));
        end
        sumrate_ZF_temp(s,i) = sum(rate_ZF(:,s,i));

    end

    % Calculate average over Monte Carlo iterations
    sumrate_LMMSE(s) = mean(sumrate_LMMSE_temp(s,:));
    sumrate_MRC(s) = mean(sumrate_MRC_temp(s,:));
    sumrate_ZF(s) = mean(sumrate_ZF_temp(s,:));
    disp(SNRdB(s));
end

% Plot simulation results
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure;
hold on; box on; grid on;
plot(SNRdB, sumrate_LMMSE, 'k--', 'LineWidth', 2);
plot(SNRdB, sumrate_MRC, 'r-o', 'LineWidth', 2);
plot(SNRdB, sumrate_ZF, 'b-', 'LineWidth', 2);
xlabel('SNR [dB]','Interpreter','latex');
ylabel('Sum rate [bit/symbol]','Interpreter','latex');
set(gca,'fontsize',16);
legend({'LMMSE','MRC','ZF'},'Interpreter','latex','Location','NorthWest');
title('Hiệu suất các kỹ thuật kết hợp thu trong Massive MIMO đa người dùng đường lên','Interpreter','latex');
grid on;
