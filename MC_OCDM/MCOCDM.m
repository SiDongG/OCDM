%% Multi-Carrier OCDM 
%% Parameter Initialization
% Use 64 sub-carriers, we spread the data symbols over a smaller set of
% subcarriers, the size of each subset has to be divisible by total
% sub-carrier length. 
K=8; %Size of subset 
N=64; %Number of Subcarrier, assume always even 
L=4; %Channel Length
Block_Num=100; %Block Number
C=4; %Len Cyclic Prefix 
M=4; %4QAM Modulation
SNR=10;
P=N+C;
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a)*(b)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);
%% Bits Generation and Symbol Generation
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Bits3=qammod(Bits2,M)*sqrt(0.5);
end
Symbols=reshape(Bits3,N,1,Block_Num);
%% Multi-Carrier Spreading 
DFnT0=zeros(K);
if mod(K,2)==0
    for m=1:K
        for n=1:K
            DFnT0(m,n)=sqrt(1/K)*exp(-1i*pi/4)*exp(1i*pi*(m-n)^2/K);
        end
    end
end
IDFnT0=DFnT0';
Symbols2=zeros(size(Symbols));
Index=1;
for a=1:N/K
    Tu=zeros(K,N);
    for b=1:K
        for c=1:N
            if c==Index
                Tu(b,c)=1;
            else
                Tu(a,c)=0;
            end
        end
        Index=Index+1;
    end
    Symbols2=Symbols2+
end




for a=1:length(u)
    for c=1:N
        if c==u(a)
            Tu(a,c)=1;
        else
            Tu(a,c)=0;
        end
    end
end
for a=1:length(b)
    for c=1:N
        if c==b(a)
            Tb(a,c)=1;
        else
            Tb(a,c)=0;
        end
    end
end