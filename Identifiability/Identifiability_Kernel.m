clear; clc; close all;
loop_Num=10;

step=4;
Start=10;
End=30;
total=zeros(1,length(Start:step:End));  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,length(Start:step:End));

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for loop=1:loop_Num
        BER=OCDM_Identifiability(SNR);            
        ratio(1,SNRdB/step-Start/step+1)=BER;
        total(1,SNRdB/step-Start/step+1)=total(1,SNRdB/step-Start/step+1)+ratio(1,SNRdB/step-Start/step+1);
    end
end

total=total/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total(:,:,1),'bx-');



set(gca,'Yscale','log');
ylim([1e-6 1]);
xlabel('SINR(dB)');
ylabel('BER');
legend('CFO,Null=1')