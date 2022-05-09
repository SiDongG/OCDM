function Bitsre=FresnelRe(M,Block_Num,N,C,~,Symbols0,H0,~)
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
%% DFnT
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
Symbols3=zeros(N,1,Block_Num);
for count=1:Block_Num
    Symbols3(:,:,count)=DFnT0*Symbols2(:,:,count);
end
%% Equalization
S=eye(N);
T=[S(2*N-P+1:N,:);S];
H=R*H0*T;
D=FFT*H*IFFT;
Symbols4=zeros(N,1,Block_Num);
for count=1:Block_Num
    Symbols4(:,:,count)=Symbols3(:,:,count)./diag(D);
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


