function BER=OCDM_Identifiability(SNR,Mode)
Block_Num=100; %Block Number
M=4; %Constellation
L=4; %Channel Order
N=16; %Block Size
K=15; %Actual Used Sub-carrier 
w=randi([-100,100])/100; %Normalized CFO
c1=1/(2*N);
c2=1/(2*N);
% Mode=1;
% SNR=10;
%% Channel Generation
h=(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
D=diag(fft(h,N));
%% Matrices
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

Df=zeros(N); 
for n=1:N
    Df(n,n)=exp(1i*2*pi*w*(n-1)/N);
end

Vc1=zeros(N);
Vc2=zeros(N);
for c=1:N
    Vc1(c,c)=exp(-1i*2*pi*c1*(c-1)^2);
    Vc2(c,c)=exp(-1i*2*pi*c2*(c-1)^2);
end

A = FFT*Vc1'*IFFT*Vc2';

Tzp = zeros(N,K);
Tzp(1:K,1:K) = eye(K);
%% Modulation
Bits=randi(0:1,[1,K*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
Bits3=qammod(Bits2,M)*sqrt(0.5);
Symbols=reshape(Bits3,K,1,Block_Num);
%% Channel and CFO Matrix
nr=randn(N,1,Block_Num);
ni=randn(N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Symbols2=zeros(N,1,Block_Num);

if Mode==0
    for a=1:Block_Num
        Symbols2(:,:,a)=Df*IFFT*D*A*Tzp*Symbols(:,:,a)+Noise(:,:,a);
    end
else
    for a=1:Block_Num
        Symbols2(:,:,a)=Df*IFFT*D*Tzp*Symbols(:,:,a)+Noise(:,:,a);
    end
end

%% Calculate Covariance Matrix 
Ryy=zeros(N);
for c=1:Block_Num
    Ryy=Ryy+Symbols2(:,:,c)*Symbols2(:,:,c)';
end
Ryy=Ryy/Block_Num;

%% CFO Synchronization
[U,V,W]=svd(IFFT*D*A*Tzp);
J=zeros(201,1);
Index=0;

if Mode==0
    for w2=-1:0.01:1
        Dff=zeros(N); 
        for n=1:N
            Dff(n,n)=exp(-1i*2*pi*w2*(n-1)/N);
        end
        Index=Index+1;
        for k=1:N-K
            U1=U';
            LNS=U1(N-k,:);
            J(Index)=J(Index)+LNS*inv(Dff)*Ryy*Dff*LNS';
        end
    end
else
    for w2=-1:0.01:1
        Dff=zeros(N); 
        for n=1:N
            Dff(n,n)=exp(-1i*2*pi*w2*(n-1)/N);
        end
        Index=Index+1;
        for k=K:N-1
            f=zeros(N,1);
            for a=1:N
                f(a)=exp(1i*(a-1)*2*pi*k/N);
            end
            J(Index)=J(Index)+f'*inv(Dff)*Ryy*Dff*f;
        end
    end
end
Index=find(J==min(J));
Est_w=1-0.01*Index;
%% CFO Compensation
Dff=zeros(N); 
for n=1:N
    Dff(n,n)=exp(-1i*2*pi*Est_w*(n-1)/N);
end
Symbols3=zeros(size(Symbols2));
for count=1:Block_Num
    Symbols3(:,:,count)=Dff*Symbols2(:,:,count);
end
%% Equalization 
Symbols_5=zeros(size(Symbols));

if Mode==0
    for count=1:Block_Num
        Symbols_5(:,:,count) = qam_sphere_decoder(IFFT*D*A*Tzp,Symbols3(:,:,count),M,Symbols(:,:,count),K);
    end
else
    for count=1:Block_Num
        Symbols_5(:,:,count) = qam_sphere_decoder(IFFT*D*Tzp,Symbols3(:,:,count),M,Symbols(:,:,count),K);
    end
end

%% Demodulation
if M==4
    Symbols6=qamdemod(Symbols_5/sqrt(1/2),M);
end
Bitsre=zeros(1,K*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:K
        dec=dec2bin(Symbols6(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
%% Error Count
BER=sum(Bitsre~=Bits)/length(Bits);