clear; clc; close all;
N=64; %Number of Subcarrier
L=4; %Channel Length
Block_Num=100; %Block Number
C=4; %Len Cyclic Prefix 
P=N+C;
K=8;
loop_Num=100;

%% 
total=zeros(1,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,6);
step=2;
Start=0;
End=20;
for dB=Start:step:End
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        [Error_rate]=MCOCDM(K,N,L,Block_Num,C,SNR);
        ratio(1,dB/step-Start/2+1)=Error_rate;
        total(1,dB/step-Start/2+1)=total(1,dB/step-Start/2+1)+ratio(1,dB/step-Start/2+1);
    end
end
total=total/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total,'bx-');
set(gca,'Yscale','log');
ylim([1e-4 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('OCDM')