%% Multi-Carrier OCDM 
%% Parameter Initialization
function [Error_rate]=MCOCDM(SNR,Mode,M)
% Use 64 sub-carriers, we spread the data symbols over a smaller set of
% subcarriers, the size of each subset has to be divisible by total
% sub-carrier length. 
K=4; %Size of subset 
N=64; %Number of Subcarrier, assume always even 
L=4; %Channel Length
Block_Num=100; %Block Number
C=4; %Len Cyclic Prefix 
% SNR=100;
% Mode=1;
% M=4; %4QAM Modulation
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
for k=1:log2(M):length(Bits)
    Bits2(1+(k-1)/log2(M))=bin2dec(num2str(Bits(k:k+log2(M)-1)));
end
if M==2
    Bits3=qammod(Bits2,M);
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
for count=1:Block_Num
    Index=1;
    for a=1:N/K
        Tu=zeros(K,N);
        for b=1:K
            for c=1:N
                if c==Index
                    Tu(b,c)=1;
                else
                    Tu(b,c)=0;
                end
            end
            Index=Index+1;
        end
        Symbols2(:,:,count)=Symbols2(:,:,count)+Tu.'*IDFnT0*Symbols(Index-K:Index-1,:,count);
    end
end
%% Channel
h=(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
H0=zeros(P); %Preallocating for speed, H0 is the P by P matrix have the (i,j)th entry h(i-j)
H1=zeros(P); %Preallocating for speed, H1 is the P by P matrix have the (i,j)th entry h(P+i-j)
a=1;
while a<P+1  %generate the channel matrces
    b=1;
    while b<P+1
        if a-b<0 || a-b>L-1
            H0(a,b)=0;
        else
            H0(a,b)=h(a-b+1);
        end
        if P+a-b<0 || P+a-b>L-1
            H1(a,b)=0;
        else
            H1(a,b)=h(P+a-b+1);
        end
        b=b+1;
    end
    a=a+1;
end
H=R*H0*T;
D=FFT*H*IFFT;
%% AWGN
D2=zeros(1,N);
if Mode==2 %Interleave
    for k=1:N/4
        D2(4*k-3:4*k)=[D(k,k),D(k+16,k+16),D(k+32,k+32),D(k+48,k+48)];
    end
    D=diag(D2);
end
nr=randn(N,1,Block_Num);
ni=randn(N,1,Block_Num);
Noise=(sqrt(2)/2)*(nr+1i*ni);
Symbols3=zeros(size(Symbols2));
for count=1:Block_Num
    Symbols3(:,:,count)=D*Symbols2(:,:,count)+(1/sqrt(SNR))*Noise(:,:,count);
end
%% Despreading 
Symbols4=zeros(size(Symbols3));
G=pinv(FFT*H*IFFT); %Composite Equalization Matrix
for count=1:Block_Num
    Index=1;
    for a=1:N/K
        Tu=zeros(K,N);
        for b=1:K
            for c=1:N
                if c==Index
                    Tu(b,c)=1;
                else
                    Tu(b,c)=0;
                end
            end
            Index=Index+1;
        end
%         Symbols4(Index-K:Index-1,:,count)=DFnT0*G(Index-K:Index-1,Index-K:Index-1)*Tu*FFT*Symbols3(:,:,count);
        Symbols3_5=DFnT0*Tu*Symbols3(:,:,count);
        Symbols4(Index-K:Index-1,:,count)=ML_detector(Symbols3_5,M,DFnT0*Tu*D*Tu.'*IDFnT0);
    end
end
%% Demodulation
if M==2
    Symbols5=qamdemod(Symbols4,M);
end
if M==4
    Symbols5=qamdemod(Symbols4/sqrt(1/2),M);
end
if M==16
    Symbols5=qamdemod(Symbols4/sqrt(1/10),M);
end
if M==64
    Symbols5=qamdemod(Symbols4/sqrt(1/42),M);
end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbols5(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
%% Error Statistics
Error_rate=sum(Bits~=Bitsre)/length(Bits);

























