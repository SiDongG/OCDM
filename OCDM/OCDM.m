
clear; clc; close all;
N=64; %Number of Subcarrier
L=4; %Channel Length
Block_Num=100; %Block Number
M=4; %Modulation QAM
C=4; %Len Cyclic Prefix 
P=68;
loop_Num=10;
Equal=1;
total=zeros(1,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,6);
for dB=0:4:20
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        [Bits,Trans_Symbols]=Transmitter(M,Block_Num,N,C);
        [H0,Symbols0]=Channel(Trans_Symbols,L,N,Block_Num,SNR);
        Bitsre=Receiver(M,Block_Num,N,C,Equal,Symbols0,H0,SNR);             
        ratio(dB/4+1)=sum(Bits~=Bitsre)/(Block_Num*N);
        total(dB/4+1)=total(dB/4+1)+ratio(dB/4+1);
    end
end
total=total/loop_Num;
figure()
semilogy(0:4:20,total);

xlabel('SNR(dB)');
ylabel('Ber');