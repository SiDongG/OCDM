clear; clc; close all;
loop_Num=500;

step=4;
Start=0;
End=40;
total=zeros(1,length(Start:step:End),2);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,length(Start:step:End),2);

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for Mode=1:2
        for loop=1:loop_Num
            [BER]=AFDMCFO(SNR,Mode);            
            ratio(1,SNRdB/step-Start/2+1,Mode)=BER;
            total(1,SNRdB/step-Start/2+1,Mode)=total(1,SNRdB/step-Start/2+1,Mode)+ratio(1,SNRdB/step-Start/2+1,Mode);
        end
    end
end

total=total/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total(:,:,1),'bx-');
plot(Start:step:End,total(:,:,2),'rx-');

set(gca,'Yscale','log');
ylim([1e-5 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('CFO-ZF','CFO-MMSE')