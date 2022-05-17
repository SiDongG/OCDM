function [Bits,Symbols2]=CFOTx(M,Block_Num,N)
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

end