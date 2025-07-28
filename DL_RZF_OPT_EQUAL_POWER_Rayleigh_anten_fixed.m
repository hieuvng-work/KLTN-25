close all;
clear;
clc;

% Select the range of SNR values
SNRdB = -5;
SNR = db2pow(SNRdB);

% Select range of the number of antennas
M = 32;
K = 8;

NTry = 1000;

gain_dB = linspace(0, 20, K);      % Pathloss 0-20 dB
pathloss = 10.^(-gain_dB/10);      % Tuyến tính


totalPower = K;
% Prepare to save simulation results



% Temporary variables for Monte Carlo iterations
sumrate_ZF_temp = zeros(1,NTry);
sumrate_ZF_opt_temp = zeros(1,NTry);

% Store all Monte Carlo results for individual users

SINR_ZF = zeros(K,NTry);
rate_ZF = zeros(K,NTry);
SINR_ZF_opt = zeros(K,NTry);
rate_ZF_opt = zeros(K,NTry);
lambdaInv_all = zeros(K,NTry);
powerCoefNorm_opt_all = zeros(K,NTry);
waterLevel_all = zeros(1,NTry);
% Monte Carlo simulation for current antenna configuration

for n = 1:NTry
        % Generate Rayleigh fading channel
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);
        % H = H .* sqrt(pathloss.*shadowing); 
        H = H .* sqrt(pathloss);
        % ZF Precoding equal power
        powerCoefNorm = totalPower/K*ones(K,1);
        % P_ZF = inv(H.' * conj(H));
        P_ZF = inv(H'*H);
        P_ZF = diag(P_ZF);     
        for k = 1:K
            SINR_ZF(k,n)= SNR*powerCoefNorm(k)/real(P_ZF(k));
            rate_ZF(k,n) = log2(1+ SINR_ZF(k,n));   
        end
        sumrate_ZF_temp(n) = sum(rate_ZF(:,n));   


        % ZF Precoding optimal power
        lambdaInv = (1/SNR)*real(P_ZF).*ones(K,1);
        [powerCoefNorm_opt,waterLevel] = functionWaterfilling(totalPower, lambdaInv);      
            
        for k = 1:K
            SINR_ZF_opt(k,n)= SNR*powerCoefNorm_opt(k)/real(P_ZF(k));
            rate_ZF_opt(k,n) = log2(1+ SINR_ZF_opt(k,n));   
        end
        sumrate_ZF_opt_temp(n) = sum(rate_ZF_opt(:,n)); 
        lambdaInv_all(:,n) = lambdaInv;
        powerCoefNorm_opt_all(:,n) = powerCoefNorm_opt;
        waterLevel_all(n) = waterLevel;
end
    
% Calculate average results for this antenna configuration
sumrate_ZF = mean(sumrate_ZF_temp);
sumrate_ZF_opt = mean(sumrate_ZF_opt_temp);  
lambdaInv_avg = mean(lambdaInv_all,2);
powerCoefNorm_opt_avg = mean(powerCoefNorm_opt_all,2);
waterLevel_avg = mean(waterLevel_all);
figure(1);

subplot(2,1,1);
bar(1:K, lambdaInv_avg,'cyan');
hold on;
plot([0.5, K+0.5], [waterLevel_avg, waterLevel_avg], 'r--', 'LineWidth', 1);
hold on;

xlabel('User','Interpreter','latex','FontSize',12);
ylabel('Value','Interpreter','latex','FontSize',12);
title('$\lambda^{-1}$ For Each User','Interpreter','latex','FontSize',12);
legend('\lambda^{-1} value','water level (\mu)','Location','best','FontSize',12);
grid on;
set(gca,'fontsize',12);
for k = 1:K
    text(k, lambdaInv_avg(k)+0.001, sprintf('%.4f', lambdaInv_avg(k)), ...
         'HorizontalAlignment', 'center','FontSize',12, 'VerticalAlignment', 'bottom');
end
% ylim([0 1.5]);
y_limits = ylim;
offset = (y_limits(2) - y_limits(1)) * 0.05;
% Thêm chú thích μ bên phải đường nét đứt
text(8.5, waterLevel_avg, sprintf('\\mu = %.4f', waterLevel_avg), ...
     'HorizontalAlignment', 'left','FontSize', 12, 'Color', 'red');


subplot(2,1,2);
bar(1:K, powerCoefNorm_opt_avg);hold on;
plot([0.5, K+0.5], [waterLevel_avg, waterLevel_avg], 'r--', 'LineWidth', 1);
xlabel('User','Interpreter','latex','FontSize',12);
ylabel('Value','Interpreter','latex','FontSize',12);
title('Optimal Power Allocation ($P_{total}= K$)','Interpreter','latex','FontSize',12);
legend('Allocated Power Level', 'water level (\mu)', 'Location', 'best', 'FontSize', 12);
grid on;

set(gca,'fontsize',12);
for k = 1:K
    text(k, powerCoefNorm_opt_avg(k)+0.001, sprintf('%.4f', powerCoefNorm_opt_avg(k)), ...
         'HorizontalAlignment', 'center','FontSize', 12, 'VerticalAlignment', 'bottom');
end

ylim(y_limits);
offset1 = (y_limits(2) - y_limits(1)) * 0.05; % 5% của range
% Thêm chú thích μ bên phải đường nét đứt
text(8.5, waterLevel_avg, sprintf('\\mu = %.4f', waterLevel_avg), ...
     'HorizontalAlignment', 'left','FontSize', 12, 'Color', 'red');
% ylim([0 1.5]);



