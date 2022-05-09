
clear; clc; close all;
N=4; %Number of Subcarrier
L=2; %Channel Length
Block_Num=2; %Block Number
M=4; %Modulation QAM
C=2; %Len Cyclic Prefix 
P=N+C;
loop_Num=10;
Equal=1;
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];

% Construct Fresnel Matrix
DFnT0=zeros(N);
DFnT1=zeros(N);
if mod(N,2)==0
    for m=1:N
        for n=1:N
            DFnT0(m,n)=sqrt(1/N)*exp(-1i*pi/4)*exp(1i*pi*(m-n)^2/N);
        end
    end
end
if mod(N,2)==1
    for m=1:N
        for n=1:N
            DFnT1(m,n)=sqrt(1/N)*exp(-1i*pi/4)*exp(1i*pi*(m-n+0.5)^2/N);
        end
    end
end        
IDFnT0=DFnT0';
IDFnT1=DFnT1';

IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);
%%
total=zeros(1,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,6);
for dB=0:4:20
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        [Bits,Symbols0]=Transmitter(M,Block_Num,N,C);
        [H0,Symbols1]=Channel(Symbols0,L,N,Block_Num,SNR);
        Bitsre=Receiver(M,Block_Num,N,C,Equal,Symbols1,H0,SNR);             
        ratio(dB/4+1)=sum(Bits~=Bitsre)/(Block_Num*N*log2(M));
        total(dB/4+1)=total(dB/4+1)+ratio(dB/4+1);
    end
end
total=total/loop_Num;
figure()
semilogy(0:4:20,total);

xlabel('SNR(dB)');
ylabel('Ber');