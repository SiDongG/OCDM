clear; clc; close all;
loop_Num=500;

step=4;
Start=0;
End=20;
total=zeros(1,length(Start:step:End),3);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,length(Start:step:End),3);
total2=zeros(1,length(Start:step:End),3);  %preallocating for Speed, SNR from 0 to 20
ratio2=zeros(1,length(Start:step:End),3);
Mod=[2,4,16];

for SNRdB=Start:step:End
    disp(SNRdB);
    SNR=10^(SNRdB/10);
    for U=1:2
        M=Mod(U);
        for loop=1:loop_Num
            [Error]=MCOCDM(SNR,1,M);            
            ratio(1,SNRdB/step-Start/2+1,U)=Error;
            total(1,SNRdB/step-Start/2+1,U)=total(1,SNRdB/step-Start/2+1,U)+ratio(1,SNRdB/step-Start/2+1,U);
            [Error_int]=MCOCDM(SNR,2,M);       
            ratio2(1,SNRdB/step-Start/2+1,U)=Error_int;
            total2(1,SNRdB/step-Start/2+1,U)=total2(1,SNRdB/step-Start/2+1,U)+ratio2(1,SNRdB/step-Start/2+1,U);
        end
    end
end

total=total/loop_Num;
total2=total2/loop_Num;
figure()
box on; hold on;
plot(Start:step:End,total(:,:,1),'bx-');
plot(Start:step:End,total(:,:,2),'rx-');
% plot(Start:step:End,total(:,:,3));
plot(Start:step:End,total2(:,:,1),'gx-');
plot(Start:step:End,total2(:,:,2));
% plot(Start:step:End,total2(:,:,3));
set(gca,'Yscale','log');
ylim([1e-5 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('OCDM 2QAM','OCDM 4QAM','OCDM 2QAM In','OCDM 4QAM In')