
close all;
clear;
clc;


SNRdB = -5;
SNR = db2pow(SNRdB);


Mall = 64:1:128;
K = 8;

NTry = 1000;

gain_dB = linspace(0, 20, K);      % Pathloss 0-20 dB
pathloss = 10.^(-gain_dB/10);      % Tuyến tính

totalPower = K;


sumrate_ZF = zeros(length(Mall),1);
sumrate_ZF_opt = zeros(length(Mall),1);

sumrate_ZF_temp = zeros(length(Mall),NTry);
sumrate_ZF_opt_temp = zeros(length(Mall),NTry);



SINR_ZF = zeros(K,length(Mall),NTry);
rate_ZF = zeros(K,length(Mall),NTry);
SINR_ZF_opt = zeros(K,length(Mall),NTry);
rate_ZF_opt = zeros(K,length(Mall),NTry);


for m = 1:length(Mall)
    M = Mall(m); 
 
    for n = 1:NTry
    
        H = (randn(M,K) + 1j*randn(M,K))/sqrt(2);
        H = H .* sqrt(pathloss);   
        powerCoefNorm = totalPower/K * ones(K,1);
        P_ZF = inv(H'*H);
        % P_ZF = inv(H.'*conj(H));
        P_ZF = diag(P_ZF);
        for k = 1:K
            SINR_ZF(k,m,n)= SNR*powerCoefNorm(k)/real(P_ZF(k));
            rate_ZF(k,m,n) = log2(1+ SINR_ZF(k,m,n));   
        end
        sumrate_ZF_temp(m,n) = sum(rate_ZF(:,m,n));    
           


        lambdaInv = (1/SNR)*real(P_ZF).*ones(K,1);
        powerCoefNorm_opt = functionWaterfilling(totalPower, lambdaInv);      
        for k = 1:K
            SINR_ZF_opt(k,m,n)= SNR*powerCoefNorm_opt(k)/real(P_ZF(k));
            rate_ZF_opt(k,m,n) = log2(1+ SINR_ZF_opt(k,m,n));   
        end
        sumrate_ZF_opt_temp(m,n) = sum(rate_ZF_opt(:,m,n)); 
    end

  
    sumrate_ZF(m) = mean(sumrate_ZF_temp(m,:));
    sumrate_ZF_opt(m) = mean(sumrate_ZF_opt_temp(m,:)); 
    
end
figure(1);
hold on; box on; grid on;
plot(Mall,sumrate_ZF,'r--','LineWidth',2);
plot(Mall,sumrate_ZF_opt,'b-','LineWidth',2);
xlabel('Number of antennas $(M)$','Interpreter','latex');
ylabel('Sum rate [bit/symbol]','Interpreter','latex');
set(gca,'fontsize',16);
legend({'ZF - Equal Power','ZF - Optimal Power Allocation'},'Interpreter','latex','Location','NorthWest');
title('Hiệu suất kỹ thuật tiền mã hóa ZF khi phân bổ công suất đều và phân bổ công suất tối ưu');
xlim([Mall(1) Mall(end)]);


