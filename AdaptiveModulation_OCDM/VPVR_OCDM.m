%% Variable Power Variable Rate Adaptive Modulation OCDM
function [Error_rate,C_Adaptation]=VPVR_OCDM(N,Block_Num,Mode,Pb,SNR)
% Variable Power Variable Rate Block-wise Bit Loading OCDM using
% Multi-carrier Spreading 
% Equal==1: ZF
% Equal==2: MMSE
% Equal==3: Sphere Decoding
% Mode==1:Consecutive
% Mode==2:Rearrange
% Mode==3:Interleaving 
% Compare rate bewteen rearrange partitioning and consecutive partitioning
% clear; clc; close all;
%% Parameter Initialization (Test)
% N=256; %Number of Subcarrier, assume always even 
% L=4; %Channel Length
% Block_Num=1; %Block Number
% C=4; %Len Cyclic Prefix 
% Equal=1;
% Mode=3;
% SNR=1;
% Pb=0.01;
%% Global Parameters
Z=10;%1/Z is the threshold searching rigidity 
C=4; %Len Cyclic Prefix 
L=4;
P=N+C;
Q=4; %Size of subset \
V=1; %Channel Variance 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
N_var=1/sqrt(SNR);
K=-1.5/(log(5*Pb)); % Above Cutoff Fade Level Constraint
Power=N; %Total Power Constraint 
Power_avg=Power/N; %1 unit of power per subcarrier 
%% Matrix Initialization 
IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

DFnT0=zeros(Q);
if mod(Q,2)==0
    for m=1:Q
        for n=1:Q
            DFnT0(m,n)=sqrt(1/Q)*exp(-1i*pi/4)*exp(1i*pi*(m-n)^2/Q);
        end
    end
end
IDFnT0=DFnT0';
%% Initialize Channel Properties 
h=V*(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
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
D=FFT*R*H0*T*IFFT;
%% Block-Wise Classification
if Mode==1
    D_img=diag(abs(D));
elseif Mode==2
    D_img=sort(diag(abs(D)));
else
    A=diag(abs(D));
    D_img=zeros(N,1);
    for i=1:N/Q
        D_img(4*i-3:4*i,1)=[A(i),A(i+N/Q),A(i+2*N/Q),A(i+3*N/Q)].';
    end
end
D_img2=Power_avg*D_img/N_var^2;
D_avg=zeros(1,length(D_img2)/Q);
for k=1:N/Q
    D_avg(k)=mean(D_img2(4*k-3:4*k));
end
% D_avg=10*log10(D_avg);
% D_avg=Power_avg*D_avg/N_var^2;
%% Sub-Carrier Classification 
y0=1; %Initialize Threshold and iteratively search for power balance
yk=y0/K;
Power_alloc=zeros(size(D_avg));
QAM_alloc=zeros(size(D_avg));
Loop=0;
Broke=0;
while abs(sum(Power_alloc)-Power/Q)>Power/Q/Z   %Within 5% Accuracy
    for count=1:N/Q
        if (0<=D_avg(count)/yk) && (D_avg(count)/yk<2)
            Power_alloc(count)=0;
            QAM_alloc(count)=0;
        elseif (2<=D_avg(count)/yk) && (D_avg(count)/yk<4)
            Power_alloc(count)=Power_avg*(1/(D_avg(count)*K));
            QAM_alloc(count)=2;
        elseif (4<=D_avg(count)/yk) && (D_avg(count)/yk<16)
            Power_alloc(count)=Power_avg*(3/(D_avg(count)*K));
            QAM_alloc(count)=4;
        elseif (16<=D_avg(count)/yk) && (D_avg(count)/yk<64)
            Power_alloc(count)=Power_avg*(15/(D_avg(count)*K));
            QAM_alloc(count)=16;
        else
            Power_alloc(count)=Power_avg*(63/(D_avg(count)*K));
            QAM_alloc(count)=64;
%         elseif (128<=D_avg(count)/yk) && (D_avg(count)/yk<256)
%             Power_alloc(count)=Power_avg*(127/(D_avg(count)*K));
%             QAM_alloc(count)=128;
%         else
%             Power_alloc(count)=Power_avg*(255/(D_avg(count)*K));
%             QAM_alloc(count)=256;
        end
    end
    if sum(Power_alloc)>Power/Q+Power/Q/Z
        yk=yk+0.0001;
    elseif sum(Power_alloc)<Power/Q-Power/Q/Z && yk>0
        yk=yk-0.0001;
    else
        yk=yk;
    end
    Loop=Loop+1;
    if Loop>100000
        Broke=1;
        break
    end
%     yk
%     sum(Power_alloc)
end
%% Expanding Power and QAM Allocation 
Power_alloc2=zeros(size(D_img));
QAM_alloc2=zeros(size(D_img));
for k=1:N/Q
    Power_alloc2(4*k-3:4*k)=Power_alloc(k);
    QAM_alloc2(4*k-3:4*k)=QAM_alloc(k);
end
% figure
% bar(D_img2)
% hold on;
% plot([0,N],[2*yk,2*yk])
% plot([0,N],[4*yk,4*yk])
% plot([0,N],[16*yk,16*yk])
% plot([0,N],[64*yk,64*yk])
% figure
% bar(Power_alloc2)
% sum(Power_alloc)
%% Bit Generation 
Total_Bits=0;
for i=1:N
    if QAM_alloc2(i)~=0
        Total_Bits=Total_Bits+log2(QAM_alloc2(i));
    end
end
%% Adaptive Modulation
Nonz=nnz(QAM_alloc2);
Symbols=zeros(N,1,Block_Num);
Bits=randi(0:1,[1,Total_Bits*Block_Num]);
Index=1;
for count=1:Block_Num
    for i=1:N
        if QAM_alloc2(i)==0
            Symbols(i,:,count)=0;
        elseif QAM_alloc2(i)==2
            B=bin2dec(num2str(Bits(Index)));
            Symbols(i,:,count)=qammod(B,2);
            Index=Index+1;
        elseif QAM_alloc2(i)==4
            B=bin2dec(num2str(Bits(Index:Index+1)));
            Symbols(i,:,count)=qammod(B,4)*sqrt(0.5);
            Index=Index+2;
        elseif QAM_alloc2(i)==16
            B=bin2dec(num2str(Bits(Index:Index+3)));
            Symbols(i,:,count)=qammod(B,16)*sqrt(1/10);
            Index=Index+4;
        else
            B=bin2dec(num2str(Bits(Index:Index+5)));
            Symbols(i,:,count)=qammod(B,64)*sqrt(1/42);
            Index=Index+6;
        end
    end
end
%% Power Adaptation 
Symbols=Power_alloc2.*Symbols;
%% Multi-Carrier Spreading
Symbols2=zeros(size(Symbols));
for count=1:Block_Num
    Index=1;
    for a=1:N/Q
        Tu=zeros(Q,N);
        for b=1:Q
            for c=1:N
                if c==Index
                    Tu(b,c)=1;
                else
                    Tu(b,c)=0;
                end
            end
            Index=Index+1;
        end
        Symbols2(:,:,count)=Symbols2(:,:,count)+Tu.'*IDFnT0*Symbols(Index-Q:Index-1,:,count);
    end
end
%% Channel
nr=randn(N,1,Block_Num);
ni=randn(N,1,Block_Num);
Noise=(sqrt(2)/2)*(nr+1i*ni);
Symbols3=zeros(size(Symbols2));
for count=1:Block_Num
    Symbols3(:,:,count)=D*Symbols2(:,:,count)+N_var*Noise(:,:,count);
end
%% Despreading and Equalization 
Symbols4=zeros(size(Symbols3));
G=pinv(D); %Composite Equalization Matrix
for count=1:Block_Num
    Index=1;
    for a=1:N/Q
        Tu=zeros(Q,N);
        for b=1:Q
            for c=1:N
                if c==Index
                    Tu(b,c)=1;
                else
                    Tu(b,c)=0;
                end
            end
            Index=Index+1;
        end
        Symbols4(Index-Q:Index-1,:,count)=DFnT0*G(Index-Q:Index-1,Index-Q:Index-1)*Tu*Symbols3(:,:,count);
    end
end
%% Preparation for calculating BER per subcarrier

%% Demodulation
Bitsre=zeros(size(Bits));
Index=0;
Symbols4=Symbols4./Power_alloc2;
Error_subcarrier=zeros(1,N);
for count=1:Block_Num
    for i=1:N
        if QAM_alloc2(i)==2
            A=qamdemod(Symbols4(i,:,count),2);
            dec=dec2bin(A,1);
            for n=1:length(dec)
                Bitsre(Index+n)=str2double(dec(n));
                if Bitsre(Index+n)~=Bits(Index+n)
                    Error_subcarrier(1,i)=Error_subcarrier(1,i)+1;
                end
            end            
            Index=Index+1;
        elseif QAM_alloc2(i)==4
            A=qamdemod(Symbols4(i,:,count)/sqrt(1/2),4);
            dec=dec2bin(A,2);
            for n=1:length(dec)
                Bitsre(Index+n)=str2double(dec(n));
                if Bitsre(Index+n)~=Bits(Index+n)
                    Error_subcarrier(1,i)=Error_subcarrier(1,i)+1;
                end
            end
            Index=Index+2;
        elseif QAM_alloc2(i)==16
            A=qamdemod(Symbols4(i,:,count)/sqrt(1/10),16);
            dec=dec2bin(A,4);
            for n=1:length(dec)
                Bitsre(Index+n)=str2double(dec(n));
                if Bitsre(Index+n)~=Bits(Index+n)
                    Error_subcarrier(1,i)=Error_subcarrier(1,i)+1;
                end
            end
            Index=Index+4;
        elseif QAM_alloc2(i)==64
            A=qamdemod(Symbols4(i,:,count)/sqrt(1/42),64);
            dec=dec2bin(A,6);
            for n=1:length(dec)
                Bitsre(Index+n)=str2double(dec(n));
                if Bitsre(Index+n)~=Bits(Index+n)
                    Error_subcarrier(1,i)=Error_subcarrier(1,i)+1;
                end
            end
            Index=Index+6;
        else
            disp('')
        end
    end
end
%% Count Error Per Subcarrier
Error_rate_subcarrier=zeros(1,N);
for i=1:N
    Error_rate_subcarrier(i)=Error_subcarrier(i)/(Block_Num*log2(QAM_alloc2(i)));
end
%% Error Counting 
Error=sum(Bitsre~=Bits);
Error_rate=Error/length(Bits);
%% Capacity Calculation 
C_Nonadaptation=log2(2);
C_Adaptation=0;
for i=1:N
    if QAM_alloc2(i)==2
        C_Adaptation=C_Adaptation+log2(2)/N;
    elseif QAM_alloc2(i)==4
        C_Adaptation=C_Adaptation+log2(4)/N;
    elseif QAM_alloc2(i)==16
        C_Adaptation=C_Adaptation+log2(16)/N;
    elseif QAM_alloc2(i)==64
        C_Adaptation=C_Adaptation+log2(64)/N;
    else
        disp('')
    end
end
Effective_SNR=Power_alloc2.*D_img2;
PowerK=sum(1/yk-1./D_img);
if Broke==1
    Error_rate=1000;
end
























