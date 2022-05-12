clear; clc; close all;
%% DFTPOFDM property 
%Parameter Definition
N=64;
N2=80;%
L=4;
C=4;
Block_Num=100;
P=N+C;
M=4;
%% Initalize Symbols, Unit Power
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Symbols=qammod(Bits2,M)*sqrt(0.5);
end
if M==16
    Symbols=qammod(Bits2,M)*sqrt(1/10);
end
if M==64
    Symbols=qammod(Bits2,M)*sqrt(1/42);
end
Symbols=reshape(Symbols,N,1,Block_Num);
%% N-point DFT
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);
Symbols0=zeros(N,1,Block_Num);
for count=1:Block_Num
    Symbols0(:,:,count)=FFT*Symbols(:,:,count);
end
%% Calculate PAPR pre-DFT-precoding
Symbols_pre=zeros(N,1,Block_Num);
for count=1:Block_Num
    Symbols_pre(:,:,count)=IFFT*Symbols(:,:,count);
end
Symbols_pre=reshape(Symbols_pre,1,N*Block_Num);
Max=abs(max(Symbols_pre));
Avg=mean(abs(Symbols_pre));
PAPR=10*log10(Max/Avg);
%% Frequency Mapping
%Use a scheme with 64 contagious subcarriers and 16 zero-paddled
%subcarriers






















