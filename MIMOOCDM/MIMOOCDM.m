%% 2 by 2 MIMO OCDM
function [Bitsre,Bits]=MIMOOCDM(Tx,Rx,L,C,M,N,Block_Num,SNR,Eq)
%% Parameter Initialization
% Standard: Complex Orthogonal STBC, 4QAM, CSIR, ZF, MMSE Equalization
% Tx=2; %Number of Transmit Antenna
% Rx=2; %Number of Receive Antenna 
% L=2;  %Channel Length
% C=2;  %CP Length
% M=4;  %4-QAM
% N=4; %Block Size
% Eq=1;
% Block_Num=Tx*2; %Number of Blocks
% SNR=1000;
P=N+C;
Var_n=1/SNR; %Noise Variance

%% Mapping 
Bits=randi(0:1,[1,N*Block_Num*log2(M)*2]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Symbols=qammod(Bits2,M)*sqrt(0.5);
end
%% DFnT and DFT matrix 
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
FFT=dftmtx(N)/sqrt(N);
IFFT=conj(FFT);
%% Reconfiguring Matrix 
Symbols=reshape(Symbols,[N,1,Tx,Block_Num]);
Symbols1=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols1(:,:,count)=[Symbols(:,:,1,count);Symbols(:,:,2,count)];
end
%% Equivalent CFR Channel Matrix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
Ch=zeros(N*Tx,N*Rx);
DD=zeros(N*Tx,N*Rx);
for i=1:Tx
    for j=1:Rx
        h=(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
        a=1;
        H0=zeros(P);
        while a<P+1  %generate the channel matrces
            b=1;
            while b<P+1
                if a-b<0 || a-b>L-1
                    H0(a,b)=0;
                else
                    H0(a,b)=h(a-b+1);
                end
                b=b+1;
            end
            a=a+1;
        end
        Ch(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=R*H0*T;
        DD(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=FFT*R*H0*T*IFFT;
    end
end
%% IDFnT
Symbols2=zeros(size(Symbols1));
if mod(N,2)==0
    for count=1:Block_Num
        Symbols2(:,:,count)=[IDFnT0*Symbols1(1:N,:,count);IDFnT0*Symbols1(N+1:2*N,:,count)];
    end
end
if mod(N,2)==1
    for count=1:Block_Num
        Symbols2(:,:,count)=[IDFnT1*Symbols1(1:N,:,count);IDFnT1*Symbols1(N+1:2*N,:,count)];
    end
end
%% Channel
nr=randn(Tx*N,1,Block_Num);
ni=randn(Tx*N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Symbols3=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols3(:,:,count)=Ch*Symbols2(:,:,count)+Noise(:,:,count);
end
%% FFT
Symbols4=zeros(size(Symbols3));
for count=1:Block_Num
    Symbols4(:,:,count)=[FFT*Symbols3(1:N,:,count);FFT*Symbols3(N+1:2*N,:,count)];
end
%% Zadoff-Chu Sequence
Y=zeros(N);
if mod(N,2)==0
    for k=1:N
        Y(N-k+1,N-k+1)=exp(-1i*pi*k^2/N);
    end
end
if mod(N,2)==1
    for k=1:N
        Y(N-k+1,N-k+1)=exp(-1i*pi*k*(k+1)/N);
    end
end
Symbols5=zeros(size(Symbols4));
for count=1:Block_Num
    Symbols5(:,:,count)=[Y*Symbols4(1:N,:,count);Y*Symbols4(N+1:2*N,:,count)];
end
%% Equalization
if Eq==1
    G=inv(DD'*DD)*DD';
else
    G=DD'*inv(DD*DD'+Var_n*eye(size(DD)));
end
Symbols6=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols6(:,:,count)=G*Symbols5(:,:,count);
end
%% IFFT
Symbols7=zeros(size(Symbols6));
for count=1:Block_Num
    Symbols7(:,:,count)=[IFFT*Symbols6(1:N,:,count);IFFT*Symbols6(N+1:2*N,:,count)];
end
%% Demod 
Bitsre=zeros(size(Bits));
Symbols8=qamdemod(Symbols7/sqrt(1/10),M);
start=1;
for count=1:Block_Num
    for k=1:Tx*N
        dec=dec2bin(Symbols8(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end