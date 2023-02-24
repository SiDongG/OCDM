clear; clc; close all;
loop_Num=10;

step=4;
Start=0;
End=20;
total=zeros(1,length(Start:step:End),2);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,length(Start:step:End),2);

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for loop=1:loop_Num
        for Mode=0:1
            BER=OCDM_Identifiability(SNR,Mode);            
            ratio(1,SNRdB/step-Start/step+1,Mode+1)=BER;
            total(1,SNRdB/step-Start/step+1,Mode+1)=total(1,SNRdB/step-Start/step+1,Mode+1)+ratio(1,SNRdB/step-Start/step+1,Mode+1);
        end
    end
end

total=total/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total(:,:,1),'bx-');
plot(Start:step:End,total(:,:,2),'gx-');



set(gca,'Yscale','log');
ylim([1e-4 1]);
xlabel('SINR(dB)');
ylabel('BER');
legend('OCDM CFO,Null=1','OFDM CFO,')