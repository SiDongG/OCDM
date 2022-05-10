function [Bits,Symbols0]=Transmitter(M,Block_Num,N,C)
P=N+C;
%% Initalize Symbols, Unit Power
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Bits3=qammod(Bits2,M)*sqrt(0.5);
end
if M==16
    Bits3=qammod(Bits2,M)*sqrt(1/10);
end
if M==64
    Bits3=qammod(Bits2,M)*sqrt(1/42);
end
%% Discrete Fresnel Transform at Transmitter
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
%% IDFnT 
Symbols=zeros(N,1,Block_Num);
Block=1;
if mod(N,2)==0
    for k=1:N:length(Bits3)
        Symbols(:,:,Block)=IDFnT0*Bits3(k:k+N-1).';
        Block=Block+1;
    end
end
if mod(N,2)==1
    for k=1:N:length(Bits3)
        Symbols(:,:,Block)=IDFnT1*Bits3(k:k+N-1).';
        Block=Block+1;
    end
end
%% Cyclic Prefix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
Symbols0=zeros(P,1,Block_Num);
for a=1:Block_Num
    Symbols0(:,:,a)=T*Symbols(:,:,a);
end
end









