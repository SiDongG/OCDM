function Bitsre=Receiver(M,Block_Num,N,C,Equal,Symbols0,H0,SNR)
%1=ZF, 2=MMSE
P=N+C;
R=[zeros(N,P-N),eye(N)];
Symbols2=zeros(N,1,Block_Num);
%% Remove Cyclic Prefix 
for a=1:Block_Num
    Symbols2(:,:,a)=R*Symbols0(:,:,a);
end
%% Construct and Apply DFT matrix 
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);
Symbols3=zeros(size(Symbols2));
% IDFT matrix normalised by 1/sqrt(N)
for count=1:Block_Num
    sample=FFT*Symbols2(:,:,count); 
    Symbols3(:,:,count)=sample;
end
%% Construct and Compensate Zadoff-Chu Sequence
Y=zeros(N);
if mod(N,2)==0
    for k=1:N
        Y(k,k)=exp(-1i*pi*k^2/N);
    end
end
if mod(N,2)==1
    for k=1:N
        Y(k,k)=exp(-1i*pi*k*(k-1)/N);
    end
end
%% Construct Equalization matrix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
H=R*H0*T;
V=FFT*H*FFT'; %%Diagonal Matrix 
if Equal==1
    G=inv(V);
end
if Equal==2
    G=conj(V)/(V^2+1/SNR);
end
%% Equalization
Symbols4=zeros(size(Symbols3));
for count=1:Block_Num
    Symbols4(:,:,count)=G*Y*Symbols3(:,:,count);
end
%% Demodulation
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
        dec=dec2bin(Symbols5(k,1,count));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
        end
    end
end
end















