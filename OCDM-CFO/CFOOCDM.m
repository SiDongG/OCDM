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

%% Matrix Initialization
% DFnT
DFnT0=zeros(NP);
DFnT1=zeros(NP);
if mod(NP,2)==0
    for m=1:NP
        for n=1:NP
            DFnT0(m,n)=sqrt(1/NP)*exp(-1i*pi/4)*exp(1i*pi*(m-n)^2/NP);
        end
    end
end
if mod(NP,2)==1
    for m=1:NP
        for n=1:NP
            DFnT1(m,n)=sqrt(1/NP)*exp(-1i*pi/4)*exp(1i*pi*(m-n+0.5)^2/NP);
        end
    end
end        
IDFnT0=DFnT0';
IDFnT1=DFnT1';
% DFT
IFFT=zeros(NP);
for a=1:NP
    for b=1:NP
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/NP);
    end 
end
IFFT=IFFT*1/sqrt(NP);
FFT=conj(IFFT);
% CFO
D=zeros(P,P);
for count=1:P
    D(count,count)=exp(1i*W*2*pi*(count-1)/NP);
end
%%












