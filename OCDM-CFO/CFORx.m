function Bitsre=CFORx(Block_Num,Equal,Symbols4,SNR)

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
for count=1:Block_Num
    Symbols7(:,:,count)=DFnT0*Symbols6(:,:,count);
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
H=fft(h0,N);
Dh=zeros(N);
for count=1:N
    Dh(N,N)=H(count);
end
%% Construct Equalization Matrix
if Equal==1
    G=inv(Dh);
end
if Equal==2
    G=Dh'/(Dh*Dh'+1/SNR);
end
%% Equalization
Symbols8=zeros(size(Symbols7));
for count=1:Block_Num
    Symbols8(:,:,count)=G*Symbols7(:,:,count);
end

end




