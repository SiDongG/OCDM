clear; clc; close all;
Tx=2; %Number of Transmit Antenna
Rx=2; %Number of Receive Antenna 
L=4;  %Channel Length
C=4;  %CP Length
M=4;  %4-QAM
N=64; %Block Size
Block_Num=Tx*2; %Number of Blocks
loop_Num=1000;

%%
total=zeros(1,9,2);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,9,2);
for dB=0:4:32
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        for Eq=1:2
            [Bitsre,Bits]=MIMOOCDM(Tx,Rx,L,C,M,N,Block_Num,SNR,Eq);
            ratio(1,dB/4+1,Eq)=sum(Bits~=Bitsre)/(length(Bits));
            total(1,dB/4+1,Eq)=total(1,dB/4+1,Eq)+ratio(1,dB/4+1,Eq);
        end
    end
end

total=total/loop_Num;
figure()
box on; hold on;
plot(0:4:32,total(:,:,1),'bx-');
plot(0:4:32,total(:,:,2),'rx-');
set(gca,'Yscale','log');
ylim([1e-6 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('MIMO-OCDM-ZF','MIMO-OCDM-MMSE')