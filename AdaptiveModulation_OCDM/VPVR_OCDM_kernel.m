clear; clc; close all;
N=256; %Number of Subcarrier
L=4; %Channel Length
Block_Num=100; %Block Number
C=4; %Len Cyclic Prefix 
P=N+C;
% N_var=1;
loop_Num=10;
Pb=1e-4;

%% 
step=2;
start=0;
End=18;
total=zeros(2,length(start:step:End),3);  %preallocating for Speed, SNR from 0 to 20
Capa=zeros(2,length(start:step:End),3);
BER=zeros(1,length(start:step:End),3);
for dB=start:step:End
    disp(['dB:' num2str(dB)]);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        disp(['loop:' num2str(loop)])
        for Mode=1:3
            for Equal=1:2
                [Error_rate,C_Adaptation]=VPVR_OCDM(N,L,Block_Num,Mode,Equal,Pb,SNR);
                Capa(Equal,dB/step-start/2+1,Mode)=C_Adaptation;
                total(Equal,dB/step-start/2+1,Mode)=total(Equal,dB/step-start/2+1,Mode)+Capa(Equal,dB/step-start/2+1,Mode);
                if Equal==1
                    BER(1,dB/step-start/2+1,Mode)=BER(1,dB/step-start/2+1,Mode)+Error_rate;
                end
            end
        end
    end
end

BER=BER/loop;
total=total/loop_Num;
figure()
box on; hold on;
plot(start:step:End,total(1,:,1),'bx-');
plot(start:step:End,total(2,:,1),'rx-');
plot(start:step:End,total(1,:,2),'gx-');
plot(start:step:End,total(2,:,2),'mx-');
plot(start:step:End,total(1,:,3),'yx-');
plot(start:step:End,total(2,:,3),'kx-');
xlabel('SNR(dB)');
ylabel('Capacity');
legend('Conse-ZF','Conse-MMSE','Smart-ZF','Smart-MMSE','Inter-ZF','Inter-MMSE')