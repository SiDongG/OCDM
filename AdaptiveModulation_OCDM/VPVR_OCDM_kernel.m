clear; clc; close all;
N=256; %Number of Subcarrier
L=4; %Channel Length
Block_Num=10; %Block Number
C=4; %Len Cyclic Prefix 
P=N+C;
V=1;
% N_var=1;
loop_Num=100;
Pb=1e-3;

%% 
step=2;
start=0;
End=18;  %preallocating for Speed, SNR from 0 to 20
Capa=zeros(1,length(start:step:End),3);
BER=zeros(1,length(start:step:End),3);
Drop=zeros(1,length(start:step:End),3);
for dB=start:step:End
    disp(['dB:' num2str(dB)]);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        disp(['loop:' num2str(loop)])
        for Mode=1:3
            [Error_rate,C_Adaptation]=VPVR_OCDM(N,Block_Num,Mode,Pb,SNR);
            if Error_rate==1000
                Drop(1,dB/step-start/2+1,Mode)=Drop(1,dB/step-start/2+1,Mode)+1;
            else
                Capa(1,dB/step-start/2+1,Mode)=Capa(1,dB/step-start/2+1,Mode)+C_Adaptation;
                BER(1,dB/step-start/2+1,Mode)=BER(1,dB/step-start/2+1,Mode)+Error_rate;
            end
        end
    end
end
for dB=start:step:End
    for Mode=1:3
        BER(1,dB/step-start/2+1,Mode)=BER(1,dB/step-start/2+1,Mode)/(loop_Num-Drop(1,dB/step-start/2+1,Mode));
        Capa(1,dB/step-start/2+1,Mode)=Capa(1,dB/step-start/2+1,Mode)/(loop_Num-Drop(1,dB/step-start/2+1,Mode));
    end
end
figure()
box on; hold on;
plot(start:step:End,Capa(1,:,1),'bx-');
plot(start:step:End,Capa(1,:,2),'gx-');
plot(start:step:End,Capa(1,:,3),'yx-');
xlabel('SNR(dB)');
ylabel('Capacity');
legend('Conse-ZF','Conse-MMSE','Smart-ZF')