function [H0,Symbols1]=Channel(Symbols0,L,N,Block_Num,SNR)
P=N+L;
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
nr=randn(P,1,Block_Num);
ni=randn(P,1,Block_Num);
Noise=(sqrt(2)/2)*(nr+1i*ni);

Symbols1=zeros(P,1,Block_Num);
for a=1:Block_Num
    Insertion1=Symbols0(:,:,a);
    if a==1
        Insertion2=zeros(P,1);
    else
        Insertion2=Symbols0(:,:,a-1);
    end
    Symbols1(:,:,a)=H0*Insertion1+H1*Insertion2+(1/sqrt(SNR))*Noise(:,:,a);
end
end