% Deterministic Channel 
h=[1/sqrt(2),1i/sqrt(2),1/sqrt(2),1i/sqrt(2),1/sqrt(2),1i/sqrt(2),1/sqrt(2),1i/sqrt(2)];
% This channel has 6 zeros on FFT grid located at position 3,5,7,11,13,15
w0=0.1;

N=16;
K=15;
Tzp = zeros(N,K);
Tzp(1:K,1:K) = eye(K);
Rss=eye(K);

c1=1/(2*N);
c2=1/(2*N);

IFFT=zeros(N);
for a=1:N
    for b=1:N
        IFFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N);
    end 
end
IFFT=IFFT*1/sqrt(N);
FFT=conj(IFFT);

Vc1=zeros(N);
Vc2=zeros(N);
for c=1:N
    Vc1(c,c)=exp(-1i*2*pi*c1*(c-1)^2);
    Vc2(c,c)=exp(-1i*2*pi*c2*(c-1)^2);
end

Df=zeros(N); 
for n=1:N
    Df(n,n)=exp(-1i*2*pi*w0*(n-1)/N);
end

D=diag(fft(h,16));
A = FFT*Vc1'*IFFT*Vc2';

% OCDM Cost Function
[U,V,W]=svd(IFFT*D*A*Tzp);

Ryy=Df*IFFT*D*A*Tzp*Rss*Tzp'*A'*D'*FFT*Df'+eye(N);
J=zeros(201,1);
Index=0;
for w=-1:0.01:1
    Dff=zeros(N); 
    for n=1:N
        Dff(n,n)=exp(-1i*2*pi*w*(n-1)/N);
    end
    Index=Index+1;
    for k=1:N-K
        U1=U';
        LNS=U1(N-k,:);
        J(Index)=J(Index)+LNS*inv(Dff)*Ryy*Dff*LNS';
    end
end

hold on

plot(-1:0.01:1,J)
xlabel('w')
ylabel('J(w)')
legend('Null=6','Null=1')

% OFDM Cost Function 

J2=zeros(201,1);
Index=0;
for w=-1:0.01:1
    Dff=zeros(N); 
    for n=1:N
        Dff(n,n)=exp(-1i*2*pi*w*(n-1)/N);
    end
    Index=Index+1;
    for k=K:N-1
        f=zeros(N,1);
        for a=1:N
            f(a)=exp(1i*(a-1)*2*pi*k/N);
        end
        J2(Index)=J2(Index)+f'*inv(Dff)*Ryy*Dff*f;
    end
end





