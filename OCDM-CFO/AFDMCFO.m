%% AFDM-CFO Analysis
% No Doppler 
function [BER]=AFDMCFO(SNR,Mode)
Block_Num=100; %Block Number
M=4; 
L=2;
N=4;
w=0.1;
% Mode=1;
% SNR=1;
P=N+L;
c1=1/(2*N);
c2=1/(2*N);
%% Channel Generation
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
%% Matrices
S=eye(N);
T_ccp=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
for a=1:L
    T_ccp(a,a+N-L)=exp(-2i*pi*c1*N*(2*(a+N-L)-2-N));
end
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

DFnT0=zeros(N);
if mod(N,2)==0
    for m=1:N
        for n=1:N
            DFnT0(m,n)=sqrt(1/N)*exp(-1i*pi/4)*exp(1i*pi*((m-1)-(n-1))^2/N);
        end
    end
end
IDFnT0=DFnT0';

Z=zeros(N);
for count=1:N
    Z(count,count)=exp(-1i*pi*(count-1)^2/N);
end

Vc1=zeros(N);
Vc2=zeros(N);
for c=1:N
    Vc1(c,c)=exp(-1i*2*pi*c1*(c-1)^2);
    Vc2(c,c)=exp(-1i*2*pi*c2*(c-1)^2);
end
%% Modulation
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
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
%% IDAFT (Inverse Discrete Affine Fourier Transform)
Symbols2=zeros(size(Symbols));
for a=1:Block_Num
    Symbols2(:,:,a)=IDFnT0*Symbols(:,:,a);
end
%% Channel and CFO Matrix
nr=randn(N,1,Block_Num);
ni=randn(N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Df=zeros(N);
for n=1:N
    Df(n,n)=exp(1i*2*pi*w*(n-1)/N);
end
delta=exp(-1i*pi*w*(N-1)/N)*FFT*Df*IFFT;
Symbols3=zeros(size(Symbols));
for count=1:Block_Num
    Symbols3(:,:,count)=Z*delta*FFT*R*H0*T_ccp*Symbols2(:,:,count)+Noise(:,:,count);
end
%% Equalization 
H=FFT*R*H0*T_ccp*IFFT;
if Mode==1
    G=pinv(H);
elseif Mode==2
    G=H'/(H*H'+1/SNR);
end
Symbols5=zeros(size(Symbols));
for count=1:Block_Num
    Symbols5(:,:,count)=IFFT*G*Symbols3(:,:,count);
end
% %% DAFT (Discrete Affine Fourier Transform)
% Symbols5=zeros(size(Symbols));
% for count=1:Block_Num
%     Symbols5(:,:,count)=Vc2*FFT*Vc1*Symbols4(:,:,count);
% end
%% Demodulation
if M==2
    Symbols6=qamdemod(Symbols5,M);
end
if M==4
    Symbols6=qamdemod(Symbols5/sqrt(1/2),M);
end
if M==16
    Symbols6=qamdemod(Symbols5/sqrt(1/10),M);
end
if M==64
    Symbols6=qamdemod(Symbols5/sqrt(1/42),M);
end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbols6(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
%% Error Count
BER=sum(Bitsre~=Bits)/length(Bits);














