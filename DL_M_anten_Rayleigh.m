close all;
clear;
clc;

% Select the range of SNR values
SNRdB = -10;
SNR = db2pow(SNRdB);

% Select range of the number of antennas
Mall = 8:1:128;
K = 8;
N_simul = 3000;


% Prepare to save simulation results
sumrate_LMMSE = zeros(length(Mall),1);
sumrate_MRT = zeros(length(Mall),1);
sumrate_ZF = zeros(length(Mall),1);
sumrate_RZF = zeros(length(Mall),1);

% Temporary variables for Monte Carlo iterations
sumrate_LMMSE_temp = zeros(length(Mall),N_simul);
sumrate_MRT_temp = zeros(length(Mall),N_simul);
sumrate_ZF_temp = zeros(length(Mall),N_simul);
sumrate_RZF_temp = zeros(length(Mall),N_simul);

% Store all Monte Carlo results for individual users
rate_LMMSE = zeros(K,length(Mall),N_simul);
SINR_LMMSE = zeros(K,length(Mall),N_simul);
rate_MRT = zeros(K,length(Mall),N_simul);
SINR_MRT = zeros(K,length(Mall),N_simul);
SINR_ZF = zeros(K,length(Mall),N_simul);
rate_ZF = zeros(K,length(Mall),N_simul);
SINR_RZF = zeros(K,length(Mall),N_simul);
rate_RZF = zeros(K,length(Mall),N_simul);


for m = 1:length(Mall)
    M = Mall(m); 
    % Monte Carlo simulation for current antenna configuration
    for n = 1:N_simul
        % Generate Rayleigh fading channel
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);

        powerCoefNorm = ones(K,1);
        % LMMSE Precoding
        P_LMMSE = (conj(H)*diag(powerCoefNorm)*H.'+ eye(M)/SNR)\conj(H);
        P_LMMSE = P_LMMSE*diag(1./vecnorm(P_LMMSE));
        innerProd_LMMSE = H.'*P_LMMSE;
        for k = 1:K
            SINR_LMMSE(k,m,n) = SNR*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2 ...
                /(SNR*abs(innerProd_LMMSE(k,:)).^2*powerCoefNorm-SNR*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2+1);
            rate_LMMSE(k,m,n) = log2(1+ SINR_LMMSE(k,m,n));
        end
        sumrate_LMMSE_temp(m,n) = sum(rate_LMMSE(:,m,n));

        powerCoefNorm = ones(K,1);
        % RZF Precoding
        P_RZF = conj(H)/(H.'*conj(H)+K/SNR*eye(K));
        P_RZF = P_RZF*diag(1./vecnorm(P_RZF));
        innerProd_RZF = H.'*P_RZF;

        for k = 1:K
            SINR_RZF(k,m,n) = SNR*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2 ...
                /(SNR*abs(innerProd_RZF(k,:)).^2*powerCoefNorm-SNR*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2+1);
            rate_RZF(k,m,n) = log2(1+ SINR_RZF(k,m,n));
        end
        sumrate_RZF_temp(m,n) = sum(rate_RZF(:,m,n));

        % ZF Precoding
        powerCoefNorm = ones(K,1);
        P_ZF = conj(H/(H'*H));
        P_ZF = P_ZF*diag(1./sqrt(sum(abs(P_ZF).^2,1)));
        innerProd_ZF = H.'*P_ZF;

        for k = 1:K
            SINR_ZF(k,m,n) = SNR*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2 ...
                   /(SNR*abs(innerProd_ZF(k,:)).^2*powerCoefNorm-SNR*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2+1);
            rate_ZF(k,m,n) = log2(1+ SINR_ZF(k,m,n));   
        end
        sumrate_ZF_temp(m,n) = sum(rate_ZF(:,m,n));   



        powerCoefNorm = ones(K,1);
        % MRT Precoding
        P_MRT = conj(H);
        P_MRT = P_MRT*diag(1./vecnorm(P_MRT));
        innerProd_MRT = H.'*P_MRT;
        for k = 1:K
            SINR_MRT(k,m,n) = SNR*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2 ...
                /(SNR*abs(innerProd_MRT(k,:)).^2*powerCoefNorm-SNR*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2+1);
            rate_MRT(k,m,n) = log2(1+ SINR_MRT(k,m,n));
        end
        sumrate_MRT_temp(m,n) = sum(rate_MRT(:,m,n));
    end

    % Calculate average results for this antenna configuration
    sumrate_LMMSE(m) = mean(sumrate_LMMSE_temp(m,:));
    sumrate_RZF(m) = mean(sumrate_RZF_temp(m,:));
    sumrate_ZF(m) = mean(sumrate_ZF_temp(m,:));
    sumrate_MRT(m) = mean(sumrate_MRT_temp(m,:));
    disp(M);
end

%fprintf('Simulation completed!\n\n');

% Plot sum rate comparison
figure(1);
hold on; box on; grid on;
plot(Mall,sumrate_RZF,'r+','LineWidth',2);
plot(Mall,sumrate_ZF,'g--','LineWidth',2);
plot(Mall,sumrate_MRT,'b-','LineWidth',2);
plot(Mall,sumrate_LMMSE,'k:','LineWidth',2);
xlabel('Number of antennas $(M)$','Interpreter','latex');
ylabel('Sum rate [bit/symbol]','Interpreter','latex');
title('Hiệu suất các kỹ thuật tiền mã hóa trong MIMO đa người dùng đường xuống','Interpreter','latex');
set(gca,'fontsize',16);
legend({'RZF','ZF','MRT','LMMSE'},'Interpreter','latex','Location','NorthWest');
xlim([Mall(1) Mall(end)])


