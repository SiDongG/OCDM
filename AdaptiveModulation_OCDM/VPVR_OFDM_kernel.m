clear; clc; close all;
N=64; %Number of Subcarrier
L=4; %Channel Length
Block_Num=100; %Block Number
C=4; %Len Cyclic Prefix 
P=N+C;
N_var=1;
loop_Num=10;

%% 
total=zeros(1,6,2);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,6,2);
step=2;
Start=10;
End=20;
for dB=Start:step:End
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        [Ratio_AM,Ratio]=VPVR_OFDM(N,L,Block_Num,C,SNR);
        ratio(1,dB/step-Start/2+1,1)=Ratio;
        total(1,dB/step-Start/2+1,1)=total(1,dB/step-Start/2+1,1)+ratio(1,dB/step-Start/2+1,1);
        ratio(1,dB/step-Start/2+1,2)=Ratio_AM;
        total(1,dB/step-Start/2+1,2)=total(1,dB/step-Start/2+1,2)+ratio(1,dB/step-Start/2+1,2);
    end
end
total=total/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total(:,:,1),'bx-');
plot(Start:step:End,total(:,:,2),'rx-');
set(gca,'Yscale','log');
ylim([1e-4 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('OFDM','OFDM-AM')