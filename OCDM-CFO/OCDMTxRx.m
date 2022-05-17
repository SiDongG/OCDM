function [Bits,Bitsre]=OCDMTxRx(SNR,Equal,N,L,Block_Num,M,C,W,K,L1,Pilot0,Pilot1,E)
%% Parameter Initialization
Pilot=length(Pilot0)+length(Pilot1); %Total Pilot length
NP=N+Pilot;   %Total OFDM frame length 
P=N+C+Pilot;
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
%% Symbols Initialization
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
%% Pilot Symbols Assignment
Symbols=reshape(Bits3,N,1,Block_Num);
Symbols0=zeros(NP,1,Block_Num);
for count=1:Block_Num
    Symbols0(:,:,count)=[Pilot0.';Symbols(:,:,count);Pilot1.'];
end
%% IDFnT
Symbols1=zeros(NP,1,Block_Num);
if mod(NP,2)==0
    for count=1:Block_Num
        Symbols1(:,:,count)=IDFnT0*Symbols0(:,:,count);
    end
end
if mod(NP,2)==1
    for count=1:Block_Num
        Symbols1(:,:,count)=IDFnT1*Symbols0(:,:,count);
    end
end
%% Cyclic Prefix 
S=eye(NP);
T=[S(2*NP-P+1:NP,:);S];
Symbols2=zeros(P,1,Block_Num);
for count=1:Block_Num
    Symbols2(:,:,count)=T*Symbols1(:,:,count);
end
%% Construct Channel Matrix
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
%% Construct Noise Matrix 
nr=randn(P,1,Block_Num);
ni=randn(P,1,Block_Num);
Noise=(sqrt(2)/2)*(nr+1i*ni);
%% AWGN Multipath Channel 
Symbols3=zeros(P,1,Block_Num);
for a=1:Block_Num
    Insertion1=Symbols2(:,:,a);
    if a==1
        Insertion2=zeros(P,1);
    else
        Insertion2=Symbols2(:,:,a-1);
    end
    Symbols3(:,:,a)=H0*Insertion1+H1*Insertion2+(1/sqrt(SNR))*Noise(:,:,a);
end
%% Channel Carrier Frequency Offset 
Symbols4=zeros(P,1,Block_Num);
for count=1:Block_Num
%     Symbols4(:,:,count)=exp(1i*W*2*pi*P*count/NP)*D*Symbols3(:,:,count);
    Symbols4(:,:,count)=D*Symbols3(:,:,count);
end
%% CFO Estimation
Sum=0;
for i=1:Block_Num
    for count=(L+1):C
        Sum=Sum+conj(Symbols4(count,1,i))*Symbols4(count+NP,1,i);
    end
end
W0=(1/(2*pi))*imag(log(Sum));
%% Remove Cyclic Prefix 
R=[zeros(NP,P-NP),eye(NP)];
Symbols5=zeros(NP,1,Block_Num);
for a=1:Block_Num
    Symbols5(:,:,a)=R*Symbols4(:,:,a);
end
%% CFO Compensation 
% CFO Compensation Matrix 
D0=zeros(NP,NP);
for count=1:NP
    D0(count,count)=exp(-1i*W0*2*pi*(count-1)/NP);
end

Symbols6=zeros(NP,1,Block_Num);
for count=1:Block_Num
    Symbols6(:,:,count)=exp(-1i*2*pi*W0*C/NP)*D0*Symbols5(:,:,count);
end
%% DFT-Zadoff-Chu Option of Converting back to Frequency Domain
% %% DFT Matrix
% Symbols7=zeros(size(Symbols6));
% % IDFT matrix normalised by 1/sqrt(N)
% for count=1:Block_Num
%     Symbols7(:,:,count)=FFT*Symbols6(:,:,count);
% end
% %% Construct and Compensate Zadoff-Chu Sequence
% Y=zeros(NP);
% if mod(NP,2)==0
%     for k=1:NP
%         Y(NP-k+1,NP-k+1)=exp(-1i*pi*k^2/NP);
%     end
% end
% if mod(NP,2)==1
%     for k=1:NP
%         Y(NP-k+1,NP-k+1)=exp(-1i*pi*k*(k+1)/NP);
%     end
% end
% Symbols8=zeros(size(Symbols7));
% for count=1:Block_Num
%     Symbols8(:,:,count)=Y*Symbols7(:,:,count);
% end
%% DFnT
Symbols7=zeros(size(Symbols6));
if mod(NP,2)==0
    for count=1:Block_Num
        Symbols7(:,:,count)=DFnT0*Symbols6(:,:,count);
    end
end
if mod(NP,2)==1
    for count=1:Block_Num
        Symbols7(:,:,count)=DFnT1*Symbols6(:,:,count);
    end
end
%% Channel Estimation
% Truncation
R1=[eye(K),zeros(K,NP-K)];
S=zeros(K,1,Block_Num);
for count=1:Block_Num
    S(:,:,count)=R1*Symbols7(:,:,count);
end
% Construct B Matrix 
B=zeros(K,K);
for count=1:K
    B(count,count)=E;
end
Var_h=0.01;
Var_n=1/SNR;
% Construct Eh Matrix 
Eh=zeros(K,K);
for count=1:K
    if count<L+1
        Eh(count,count)=1/Var_h;
    else
        Eh(count,count)=0;
    end
end
% Channel Estimation
h0=1/Var_n*inv(Eh+1/Var_n*(B'*B))*B'*S(:,:,1);
h0=h0(1:L);
%% Construct Diagonal CFR Matrix
a=1;
H=zeros(P);
while a<P+1 
    b=1;
    while b<P+1
        if a-b<0 || a-b>L-1
            H(a,b)=0;
        else
            H(a,b)=h0(a-b+1);
        end
        b=b+1;
    end
    a=a+1;
end
HH=R*H*T;
Dh=FFT*HH*IFFT;
%% Construct Equalization Matrix
if Equal==1
    G=inv(HH);
end
if Equal==2
    G=HH'/(HH*HH'+1/SNR);
end
%% Equalization
Symbols8=zeros(size(Symbols7));
for count=1:Block_Num
    Symbols8(:,:,count)=G*Symbols7(:,:,count);
end
% %% Apply IDFT if DFT Zadoff-Chu Option used
% Symbols9=zeros(size(Symbols8));
% for count=1:Block_Num
%     Symbols9(:,:,count)=IFFT*Symbols8(:,:,count);
% end
%% Pilot Removal 
Symbols9=zeros(N,1,Block_Num);
for count=1:Block_Num
    Symbols9(:,:,count)=Symbols8(K+1:NP-L1,:,count);
end
%% Demodulation
if M==4
    Symbols10=qamdemod(Symbols9/sqrt(1/2),M);
end
if M==16
    Symbols10=qamdemod(Symbols9/sqrt(1/10),M);
end
if M==64
    Symbols10=qamdemod(Symbols9/sqrt(1/42),M);
end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbols10(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
end