close all;
clear;
clc;

% Select the range of SNR values
SNRdB = -10:1:30;
SNR = db2pow(SNRdB);

% Select range of the number of antennas
M = 8;
K = 8;
N_simul = 3000;

% gain_dB = linspace(0, 20, K);      % Pathloss 0-20 dB
% pathloss = 10.^(-gain_dB/10);      % Tuyến tính

% Prepare to save simulation results
sumrate_LMMSE = zeros(length(SNR),1);
sumrate_MRT = zeros(length(SNR),1);
sumrate_ZF = zeros(length(SNR),1);
sumrate_RZF = zeros(length(SNR),1);

% Temporary variables for Monte Carlo iterations
sumrate_LMMSE_temp = zeros(length(SNR),N_simul);
sumrate_MRT_temp = zeros(length(SNR),N_simul);
sumrate_ZF_temp = zeros(length(SNR),N_simul);
sumrate_RZF_temp = zeros(length(SNR),N_simul);

% Store all Monte Carlo results for individual users
rate_LMMSE = zeros(K,length(SNR),N_simul);
SINR_LMMSE = zeros(K,length(SNR),N_simul);
rate_MRT = zeros(K,length(SNR),N_simul);
SINR_MRT = zeros(K,length(SNR),N_simul);
SINR_ZF = zeros(K,length(SNR),N_simul);
rate_ZF = zeros(K,length(SNR),N_simul);
SINR_RZF = zeros(K,length(SNR),N_simul);
rate_RZF = zeros(K,length(SNR),N_simul);


for s = 1:length(SNRdB)
    % Monte Carlo simulation for current antenna configuration
    for n = 1:N_simul
        % Generate Rayleigh fading channel
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);

        powerCoefNorm = ones(K,1);
        % LMMSE Precoding
        P_LMMSE = (conj(H)*diag(powerCoefNorm)*H.'+ eye(M)/SNR(s))\conj(H);
        P_LMMSE = P_LMMSE*diag(1./vecnorm(P_LMMSE));
        innerProd_LMMSE = H.'*P_LMMSE;
        for k = 1:K
            SINR_LMMSE(k,s,n) = SNR(s)*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2 ...
                /(SNR(s)*abs(innerProd_LMMSE(k,:)).^2*powerCoefNorm-SNR(s)*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2+1);
            rate_LMMSE(k,s,n) = log2(1+ SINR_LMMSE(k,s,n));
        end
        sumrate_LMMSE_temp(s,n) = sum(rate_LMMSE(:,s,n));

        powerCoefNorm = ones(K,1);
        % RZF Precoding
        P_RZF = conj(H)/(H.'*conj(H)+K/SNR(s)*eye(K));
        P_RZF = P_RZF*diag(1./vecnorm(P_RZF));
        innerProd_RZF = H.'*P_RZF;

        for k = 1:K
            SINR_RZF(k,s,n) = SNR(s)*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2 ...
                /(SNR(s)*abs(innerProd_RZF(k,:)).^2*powerCoefNorm-SNR(s)*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2+1);
            rate_RZF(k,s,n) = log2(1+ SINR_RZF(k,s,n));
        end
        sumrate_RZF_temp(s,n) = sum(rate_RZF(:,s,n));

        % ZF Precoding
        powerCoefNorm = ones(K,1);
        P_ZF = conj(H/(H'*H));
        P_ZF = P_ZF*diag(1./sqrt(sum(abs(P_ZF).^2,1)));
        innerProd_ZF = H.'*P_ZF;

        for k = 1:K
            SINR_ZF(k,s,n) = SNR(s)*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2 ...
                   /(SNR(s)*abs(innerProd_ZF(k,:)).^2*powerCoefNorm-SNR(s)*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2+1);
            rate_ZF(k,s,n) = log2(1+ SINR_ZF(k,s,n));   
        end
        sumrate_ZF_temp(s,n) = sum(rate_ZF(:,s,n));   



        powerCoefNorm = ones(K,1);
        % MRT Precoding
        P_MRT = conj(H);
        P_MRT = P_MRT*diag(1./vecnorm(P_MRT));
        innerProd_MRT = H.'*P_MRT;
        for k = 1:K
            SINR_MRT(k,s,n) = SNR(s)*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2 ...
                /(SNR(s)*abs(innerProd_MRT(k,:)).^2*powerCoefNorm-SNR(s)*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2+1);
            rate_MRT(k,s,n) = log2(1+ SINR_MRT(k,s,n));
        end
        sumrate_MRT_temp(s,n) = sum(rate_MRT(:,s,n));
    end

    % Calculate average results for this antenna configuration
    sumrate_LMMSE(s) = mean(sumrate_LMMSE_temp(s,:));
    sumrate_RZF(s) = mean(sumrate_RZF_temp(s,:));
    sumrate_ZF(s) = mean(sumrate_ZF_temp(s,:));
    sumrate_MRT(s) = mean(sumrate_MRT_temp(s,:));
    disp(SNRdB(s));
end

%fprintf('Simulation completed!\n\n');

% Plot sum rate comparison
figure(1);
hold on; box on; grid on;
plot(SNRdB,sumrate_RZF,'r+','LineWidth',2);
plot(SNRdB,sumrate_ZF,'g--','LineWidth',2);
plot(SNRdB,sumrate_MRT,'b-','LineWidth',2);
plot(SNRdB,sumrate_LMMSE,'k:','LineWidth',2);
xlabel('SNR [dB]','Interpreter','latex');
ylabel('Sum rate [bit/symbol]','Interpreter','latex');
title('Hiệu suất các kỹ thuật tiền mã hóa trong Massive MIMO đa người dùng đường xuống','Interpreter','latex');
set(gca,'fontsize',16);
legend({'RZF','ZF','MRT','LMMSE'},'Interpreter','latex','Location','NorthWest');
xlim([SNRdB(1) SNRdB(end)])

% %% MIMO Precoding Simulation: Sum Rate & Outage per User
% close all; clear; clc;
% 
% % Simulation parameters
% SNRdB = 10;
% SNR = db2pow(SNRdB);
% Mall = 64:1:128;
% K = 16;
% N_simul = 2000;
% gamma_th_dB = 10;               % SINR threshold in dB
% gamma_th = 10^(gamma_th_dB/10);% SINR threshold in linear
% sigma_e = 0.5;
% % Preallocate result arrays
% sumrate_MRT = zeros(length(Mall),1);
% sumrate_ZF = zeros(length(Mall),1);
% sumrate_RZF = zeros(length(Mall),1);
% sumrate_LMMSE = zeros(length(Mall),1);
% 
% outage_user_MRT = zeros(K,length(Mall));
% outage_user_ZF = zeros(K,length(Mall));
% outage_user_RZF = zeros(K,length(Mall));
% outage_user_LMMSE = zeros(K,length(Mall));
% 
% for m = 1:length(Mall)
%     M = Mall(m);
% 
%     % Temp arrays
%     rate_temp_MRT = zeros(K, N_simul);
%     rate_temp_ZF = zeros(K, N_simul);
%     rate_temp_RZF = zeros(K, N_simul);
%     rate_temp_LMMSE = zeros(K, N_simul);
% 
%     outage_count_MRT = zeros(K,1);
%     outage_count_ZF = zeros(K,1);
%     outage_count_RZF = zeros(K,1);
%     outage_count_LMMSE = zeros(K,1);
% 
%     for n = 1:N_simul
%         H = (randn(M, K) + 1i*randn(M, K)) * sqrt(0.5); % Rayleigh fading
%         N_error = (randn(M, K) + 1i*randn(M, K)) * sqrt(0.5); % Lỗi ước lượng
%         H_est = sqrt(1-sigma_e^2)*H + sigma_e*N_error;
% 
%         powerCoefNorm = ones(K,1);
% 
%         % LMMSE Precoding
%         P_LMMSE = (conj(H_est)*diag(powerCoefNorm)*H_est.'+ eye(M)/SNR) \ conj(H);
%         P_LMMSE = P_LMMSE * diag(1./vecnorm(P_LMMSE));
%         inner_LMMSE = H.' * P_LMMSE;
%         for k = 1:K
%             num = SNR * powerCoefNorm(k) * abs(inner_LMMSE(k,k))^2;
%             den = SNR * (abs(inner_LMMSE(k,:)).^2 * powerCoefNorm) - num + 1;
%             SINR = num / den;
%             rate_temp_LMMSE(k,n) = log2(1 + SINR);
%             outage_count_LMMSE(k) = outage_count_LMMSE(k) + (SINR < gamma_th);
%         end
% 
%         % RZF Precoding
%         P_RZF = conj(H_est) / (H_est.'*conj(H_est) + K/SNR * eye(K));
%         P_RZF = P_RZF * diag(1./vecnorm(P_RZF));
%         inner_RZF = H.' * P_RZF;
%         for k = 1:K
%             num = SNR * powerCoefNorm(k) * abs(inner_RZF(k,k))^2;
%             den = SNR * (abs(inner_RZF(k,:)).^2 * powerCoefNorm) - num + 1;
%             SINR = num / den;
%             rate_temp_RZF(k,n) = log2(1 + SINR);
%             outage_count_RZF(k) = outage_count_RZF(k) + (SINR < gamma_th);
%         end
% 
%         % ZF Precoding
%         P_ZF = conj(H_est / (H_est'*H_est));
%         P_ZF = P_ZF * diag(1./sqrt(sum(abs(P_ZF).^2,1)));
%         inner_ZF = H.' * P_ZF;
%         for k = 1:K
%             num = SNR * powerCoefNorm(k) * abs(inner_ZF(k,k))^2;
%             den = SNR * (abs(inner_ZF(k,:)).^2 * powerCoefNorm) - num + 1;
%             SINR = num / den;
%             rate_temp_ZF(k,n) = log2(1 + SINR);
%             outage_count_ZF(k) = outage_count_ZF(k) + (SINR < gamma_th);
%         end
% 
%         % MRT Precoding
%         P_MRT = conj(H_est);
%         P_MRT = P_MRT * diag(1./vecnorm(P_MRT));
%         inner_MRT = H.' * P_MRT;
%         for k = 1:K
%             num = SNR * powerCoefNorm(k) * abs(inner_MRT(k,k))^2;
%             den = SNR * (abs(inner_MRT(k,:)).^2 * powerCoefNorm) - num + 1;
%             SINR = num / den;
%             rate_temp_MRT(k,n) = log2(1 + SINR);
%             outage_count_MRT(k) = outage_count_MRT(k) + (SINR < gamma_th);
%         end
%     end
% 
%     % Average Sum Rate
%     sumrate_MRT(m) = mean(sum(rate_temp_MRT, 1));
%     sumrate_ZF(m) = mean(sum(rate_temp_ZF, 1));
%     sumrate_RZF(m) = mean(sum(rate_temp_RZF, 1));
%     sumrate_LMMSE(m) = mean(sum(rate_temp_LMMSE, 1));
% 
%     % Per-user Outage Probabilities
%     outage_user_MRT(:,m) = outage_count_MRT / N_simul;
%     outage_user_ZF(:,m) = outage_count_ZF / N_simul;
%     outage_user_RZF(:,m) = outage_count_RZF / N_simul;
%     outage_user_LMMSE(:,m) = outage_count_LMMSE / N_simul;
% 
%     fprintf("Done M = %d\n", M);
% end
% 
% %% Plot Sum Rate
% figure;
% hold on; grid on; box on;
% plot(Mall, sumrate_RZF, 'r-o', 'LineWidth', 1);
% plot(Mall, sumrate_ZF, 'g-v', 'LineWidth', 1);
% plot(Mall, sumrate_MRT, 'b--', 'LineWidth', 1);
% plot(Mall, sumrate_LMMSE, 'k-', 'LineWidth', 1);
% xlabel('Number of BS antennas M');
% ylabel('Sum Rate [bps/Hz]');
% title('Sum Rate vs. Antennas');
% legend({'RZF','ZF','MRT','LMMSE'}, 'Location','NorthWest');
% 
% %% Plot Outage Probability (example: User 1)
% figure;
% hold on; grid on; box on;
% plot(Mall, outage_user_RZF(1,:), 'r-o', 'LineWidth', 1); hold on;
% plot(Mall, outage_user_ZF(1,:), 'c-v', 'LineWidth', 1);
% plot(Mall, outage_user_MRT(1,:), 'b--', 'LineWidth', 1);
% plot(Mall, outage_user_LMMSE(1,:), 'k-', 'LineWidth', 1);
% xlabel('Number of BS antennas M');
% title('Outage Probability')
% legend({'RZF','ZF','MRT','LMMSE'}, 'Location','NorthEast');

% % Performance comparison at specific antenna numbers
% figure(2);
% selected_M = [8,64];
% methods = {'LMMSE', 'RZF', 'ZF', 'MRT'};
% colors = {'k', 'r', 'g', 'b'};
% 
% for i = 1:length(selected_M)
%     M_idx = find(Mall == selected_M(i));
%     if ~isempty(M_idx)
%         subplot(2,2,i);
%         sumrates = [sumrate_LMMSE(M_idx), sumrate_RZF(M_idx), sumrate_ZF(M_idx), sumrate_MRT(M_idx)];
%         bar(sumrates, 'FaceColor', 'flat');
%         for j = 1:4
%             text(j, sumrates(j)+0.00001, sprintf('%.2f', sumrates(j)), ...
%                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%         end
%         set(gca, 'XTickLabel', methods);
%         ylabel('Sum rate [bit/symbol]','Interpreter','latex');
%         title(sprintf('M = %d antennas', selected_M(i)),'Interpreter','latex');
%         grid on;
%         set(gca,'fontsize',12);
%     end
% end

% % Additional analysis plots
% figure(2);
% subplot(2,2,1);
% hold on; box on; grid on;
% plot(Mall, rate_LMMSE');
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('Rate per user [bit/symbol]','Interpreter','latex');
% title('LMMSE - Rate per User','Interpreter','latex');
% set(gca,'fontsize',12);
% legend(arrayfun(@(k) sprintf('User %d', k), 1:K, 'UniformOutput', false), 'Location', 'best');
% 
% subplot(2,2,2);
% hold on; box on; grid on;
% plot(Mall, rate_RZF');
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('Rate per user [bit/symbol]','Interpreter','latex');
% title('RZF - Rate per User','Interpreter','latex');
% set(gca,'fontsize',12);
% legend(arrayfun(@(k) sprintf('User %d', k), 1:K, 'UniformOutput', false), 'Location', 'best');
% 
% subplot(2,2,3);
% hold on; box on; grid on;
% plot(Mall, rate_ZF');
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('Rate per user [bit/symbol]','Interpreter','latex');
% title('ZF - Rate per User','Interpreter','latex');
% set(gca,'fontsize',12);
% legend(arrayfun(@(k) sprintf('User %d', k), 1:K, 'UniformOutput', false), 'Location', 'best');
% 
% subplot(2,2,4);
% hold on; box on; grid on;
% plot(Mall, rate_MRT');
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('Rate per user [bit/symbol]','Interpreter','latex');
% title('MRT - Rate per User','Interpreter','latex');
% set(gca,'fontsize',12);
% legend(arrayfun(@(k) sprintf('User %d', k), 1:K, 'UniformOutput', false), 'Location', 'best');



% % SINR comparison
% figure(4);
% subplot(2,2,1);
% plot(Mall, 10*log10(SINR_LMMSE'));
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('SINR [dB]','Interpreter','latex');
% title('LMMSE - SINR per User','Interpreter','latex');
% grid on;
% set(gca,'fontsize',12);
% 
% subplot(2,2,2);
% plot(Mall, 10*log10(SINR_RZF'));
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('SINR [dB]','Interpreter','latex');
% title('RZF - SINR per User','Interpreter','latex');
% grid on;
% set(gca,'fontsize',12);
% 
% subplot(2,2,3);
% plot(Mall, 10*log10(SINR_ZF'));
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('SINR [dB]','Interpreter','latex');
% title('ZF - SINR per User','Interpreter','latex');
% grid on;
% set(gca,'fontsize',12);
% 
% subplot(2,2,4);
% plot(Mall, 10*log10(SINR_MRT'));
% xlabel('Number of antennas $(M)$','Interpreter','latex');
% ylabel('SINR [dB]','Interpreter','latex');
% title('MRT - SINR per User','Interpreter','latex');
% grid on;
% set(gca,'fontsize',12);




% close all;
% clear;
% clc;
% 
% %% Thêm các tham số mới cho xác suất dừng
% SNRdB = -10;                     % SNR ở dạng dB
% SNR = db2pow(SNRdB);             % Chuyển sang tỉ lệ tuyến tính
% gamma_th_dB = 0;                 % Ngưỡng SINR cho xác suất dừng (dB)
% gamma_th = db2pow(gamma_th_dB);  % Ngưỡng SINR dạng tuyến tính
% 
% Mall = 8:1:64;                 % Dải số anten trạm gốc
% K = 8;                           % Số người dùng
% N_simul = 10000;                  % Số lần mô phỏng Monte Carlo
% 
% %% Khởi tạo biến cho xác suất dừng
% outage_prob_LMMSE = zeros(length(Mall),1);
% outage_prob_MRT = zeros(length(Mall),1);
% outage_prob_ZF = zeros(length(Mall),1);
% outage_prob_RZF = zeros(length(Mall),1);
% 
% for m = 1:length(Mall)
%     M = Mall(m); 
% 
%     % Biến tạm đếm sự kiện dừng
%     outage_count_LMMSE = 0;
%     outage_count_MRT = 0;
%     outage_count_ZF = 0;
%     outage_count_RZF = 0;
% 
%     for n = 1:N_simul
%         % Tạo kênh truyền Rayleigh
%         H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);
%         powerCoefNorm = ones(K,1);
% 
%         %% LMMSE Precoding
%         P_LMMSE = (conj(H)*diag(powerCoefNorm)*H.'+ eye(M)/SNR)\conj(H);
%         P_LMMSE = P_LMMSE*diag(1./vecnorm(P_LMMSE));
%         innerProd_LMMSE = H.'*P_LMMSE;
% 
%         % Kiểm tra xác suất dừng cho từng người dùng
%         user_outage_LMMSE = false(K,1);
%         for k = 1:K
%             SINR_k = SNR*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2 ...
%                 /(SNR*sum(abs(innerProd_LMMSE(k,:)).^2.*powerCoefNorm.' - SNR*powerCoefNorm(k)*abs(innerProd_LMMSE(k,k))^2 + 1));
%             user_outage_LMMSE(k) = (SINR_k < gamma_th);
%         end
%         if any(user_outage_LMMSE)
%             outage_count_LMMSE = outage_count_LMMSE + 1;
%         end
% 
%         %% RZF Precoding
%         P_RZF = conj(H)/(H.'*conj(H)+K/SNR*eye(K));
%         P_RZF = P_RZF*diag(1./vecnorm(P_RZF));
%         innerProd_RZF = H.'*P_RZF;
% 
%         % Kiểm tra xác suất dừng
%         user_outage_RZF = false(K,1);
%         for k = 1:K
%             SINR_k = SNR*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2 ...
%                 /(SNR*sum(abs(innerProd_RZF(k,:)).^2.*powerCoefNorm.' - SNR*powerCoefNorm(k)*abs(innerProd_RZF(k,k))^2 + 1));
%             user_outage_RZF(k) = (SINR_k < gamma_th);
%         end
%         if any(user_outage_RZF)
%             outage_count_RZF = outage_count_RZF + 1;
%         end
% 
%         %% ZF Precoding
%         P_ZF = conj(H/(H'*H));
%         P_ZF = P_ZF*diag(1./sqrt(sum(abs(P_ZF).^2,1)));
%         innerProd_ZF = H.'*P_ZF;
% 
%         % Kiểm tra xác suất dừng
%         user_outage_ZF = false(K,1);
%         for k = 1:K
%             SINR_k = SNR*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2 ...
%                    /(SNR*sum(abs(innerProd_ZF(k,:)).^2.*powerCoefNorm.' - SNR*powerCoefNorm(k)*abs(innerProd_ZF(k,k))^2 + 1));
%             user_outage_ZF(k) = (SINR_k < gamma_th);
%         end
%         if any(user_outage_ZF)
%             outage_count_ZF = outage_count_ZF + 1;
%         end
% 
%         %% MRT Precoding
%         P_MRT = conj(H);
%         P_MRT = P_MRT*diag(1./vecnorm(P_MRT));
%         innerProd_MRT = H.'*P_MRT;
% 
%         % Kiểm tra xác suất dừng
%         user_outage_MRT = false(K,1);
%         for k = 1:K
%             SINR_k = SNR*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2 ...
%                 /(SNR*sum(abs(innerProd_MRT(k,:)).^2.*powerCoefNorm.' - SNR*powerCoefNorm(k)*abs(innerProd_MRT(k,k))^2 + 1));
%             user_outage_MRT(k) = (SINR_k < gamma_th);
%         end
%         if any(user_outage_MRT)
%             outage_count_MRT = outage_count_MRT + 1;
%         end
%     end
% 
%     % Tính xác suất dừng cho từng kỹ thuật
%     outage_prob_LMMSE(m) = outage_count_LMMSE / N_simul;
%     outage_prob_RZF(m) = outage_count_RZF / N_simul;
%     outage_prob_ZF(m) = outage_count_ZF / N_simul;
%     outage_prob_MRT(m) = outage_count_MRT / N_simul;
% 
%     fprintf('M = %d, Outage Prob: LMMSE=%.4f, RZF=%.4f, ZF=%.4f, MRT=%.4f\n', ...
%         M, outage_prob_LMMSE(m), outage_prob_RZF(m), outage_prob_ZF(m), outage_prob_MRT(m));
% end
% 
% %% Vẽ đồ thị xác suất dừng
% figure;
% hold on; box on; grid on;
% semilogy(Mall, outage_prob_RZF, 'r', 'LineWidth', 2);
% semilogy(Mall, outage_prob_ZF, 'g', 'LineWidth', 2);
% semilogy(Mall, outage_prob_MRT, 'b', 'LineWidth', 2);
% semilogy(Mall, outage_prob_LMMSE, 'k:', 'LineWidth', 2);
% xlabel('Number of antennas (M)');
% ylabel('Outage Probability');
% title('Xác suất dừng của các kỹ thuật tiền mã hóa');
% set(gca, 'YScale', 'log');
% legend({'RZF', 'ZF', 'MRT', 'LMMSE'}, 'Location', 'NorthEast');
% 
