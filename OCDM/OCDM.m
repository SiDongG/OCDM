
clear; clc; close all;
N=64; %Number of Subcarrier
L=4; %Channel Length
Block_Num=100; %Block Number
%M=4; %Modulation QAM
C=4; %Len Cyclic Prefix 
P=N+C;
loop_Num=1000;
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];

%%
total=zeros(1,11,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,11,6);
for dB=0:4:40
    disp(dB);
    SNR=10^(dB/10);
    count=1;
    for Equal=1:2
        for U=1:3
            M=4^U;
            for loop=1:loop_Num
                [Bits,Symbols0]=Transmitter(M,Block_Num,N,C);
                [H0,Symbols1]=Channel(Symbols0,L,N,Block_Num,SNR);
                Bitsre=Receiver(M,Block_Num,N,C,Equal,Symbols1,H0,SNR);             
                ratio(1,dB/4+1,count)=sum(Bits~=Bitsre)/(Block_Num*N*log2(M));
                total(1,dB/4+1,count)=total(1,dB/4+1,count)+ratio(1,dB/4+1,count);
            end
            count=count+1;
        end
    end
end
total=total/loop_Num;
figure()
box on; hold on;
plot(0:4:40,total(:,:,1),'bx-');
plot(0:4:40,total(:,:,2),'rx-');
plot(0:4:40,total(:,:,3),'gx-');
plot(0:4:40,total(:,:,4),'b-');
plot(0:4:40,total(:,:,5),'r-');
plot(0:4:40,total(:,:,6),'g-');
set(gca,'Yscale','log');
ylim([1e-4 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('OCDM/ZF 4QAM','OCDM/ZF 16QAM','OCDM/ZF 64QAM','OCDM/MMSE 4QAM','OCDM/MMSE 16QAM','OCDM/MMSE 64QAM')