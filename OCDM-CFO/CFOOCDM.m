%% Channel Estimation for OCDM w CFO 
% Author: Sidong Guo
% Date: May 19th 2022
%% Parameter Initialization
clear; clc; close all;
N=4;   %Number of Subcarrier
L=4;    %Channel Length
Block_Num=100; %Block Number
M=4;   %Modulation QAM
C=8;    %Len Cyclic Prefix, for the CFO estimation scheme, C>L
W=0.3;  %Subcarrier Frequency Offset 
SNR=1e7;
E=3;  %Total Energy for Pilots 
Pilot0=[E,0,0,0,0,0,0,0];
Pilot1=[0,0,0,0,0,0,0,0];
Equal=1;
K=length(Pilot0);
L1=length(Pilot1);
Pilot=length(Pilot0)+length(Pilot1); %Total Pilot length
NP=N+Pilot;   %Total OFDM frame length 
P=N+C+Pilot;
loop_Num=100;

%% Simulation
total=zeros(1,11,3);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,11,3);
for dB=0:2:20
    disp(dB);
    SNR=10^(dB/10);
    count=1;
    for U=1:3
        M=4^U;
        for loop=1:loop_Num
            [Bits,Bitsre]=OCDMTxRx(SNR,Equal,N,L,Block_Num,M,C,W,K,L1,Pilot0,Pilot1,E);     
            ratio(1,dB/2+1,count)=sum(Bits~=Bitsre)/(Block_Num*N*log2(M));
            total(1,dB/2+1,count)=total(1,dB/2+1,count)+ratio(1,dB/2+1,count);
        end
        count=count+1;
    end
end
total=total/loop_Num;
figure()
box on; hold on;
plot(0:2:20,total(:,:,1),'bx-');
plot(0:2:20,total(:,:,2),'rx-');
plot(0:2:20,total(:,:,3),'gx-');

set(gca,'Yscale','log');
ylim([1e-4 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('OCDM/ZF 4QAM','OCDM/ZF 16QAM','OCDM/ZF 64QAM')












