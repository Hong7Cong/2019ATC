clear all;  clc;
NN= [20:40:420];
MM = 6:3:9;
T = 100;
iter = 1;
betaE = 1;
AN_ratio = 0;
Qtotal = 10^1;%10dB
Q0max = AN_ratio*Qtotal;
Qmax = (1-AN_ratio)*Qtotal;
Pmax = 10^0.1;%1dB
RS_sum_cls1 = zeros(length(NN),length(MM));
RS_1_cls1 = zeros(length(NN),length(MM));
RS_2_cls1 = zeros(length(NN),length(MM));
RS_sum_cls1_aprx = zeros(length(NN),length(MM));
RS_1_cls1_aprx = zeros(length(NN),length(MM));
RS_2_cls1_aprx = zeros(length(NN),length(MM));
for i=1:length(MM)
    fprintf('%d of %d\n',i,length(MM));
    PdB = zeros(MM(i),2);
    QdB = zeros(MM(i),2);
    Q0dB = zeros(MM(i),2);
    beta = zeros(MM(i),2);
    for m=1:MM(i)
        PdB(m,1) = Pmax;
        PdB(m,2) = Pmax;
        QdB(m,1) = Qmax/MM(i)/2;
        QdB(m,2) = Qmax/MM(i)/2;
        Q0dB(m) = Q0max/MM(i);
        beta(m,1) = 1;
        beta(m,2) = 1;
    end
    tau = MM(i);
    [RSsum,RSsum_aprx,RS1,RS2,RS1_aprx,RS2_aprx,R1,...
        R11,R12,R1_aprx,R11_aprx,R12_aprx,...
        R2,R2_aprx,RE1,RE2,RE1_aprx,RE2_aprx,...
        Theta11,I111,I112,I113,I113_sim,I114,I114_sim...
        Theta12,I121,I122,I123,I123_sim,I124,I124_sim...
        H1,H2,Hh,W,z,...
        Q,P,Q0,rho,tmp1,tmp2,tmp3,tmp4,tmp5]...
= Sim_VarN_ZF_Chuyen(QdB,PdB,Q0dB,beta,betaE,tau,T,MM(i),NN,iter);

RS_sum_cls1(:,i) = RSsum(:,1);
RS_1_cls1(:,i) = RS1(:,1);
RS_2_cls1(:,i) = RS2(:,1);
RS_sum_cls1_aprx(:,i) = RSsum_aprx(:,1);
RS_1_cls1_aprx(:,i) = RS1_aprx(:,1);
RS_2_cls1_aprx(:,i) = RS2_aprx(:,1);
end
save('VarN_01_MRT',...
    'NN','MM','T','betaE',...
    'AN_ratio','Qtotal','Pmax',...
    'RS_sum_cls1','RS_1_cls1','RS_2_cls1',...
    'RS_sum_cls1_aprx','RS_1_cls1_aprx','RS_2_cls1_aprx');
fprintf('DONE!');
figure
plot(NN,RS_sum_cls1(:,1),'Displayname','M=3');
hold on;
plot(NN,RS_sum_cls1(:,2),'Displayname','M=6');
plot(NN,RS_sum_cls1(:,3),'Displayname','M=9');
scatter(NN,RS_sum_cls1_aprx(:,1));
scatter(NN,RS_sum_cls1_aprx(:,2));
scatter(NN,RS_sum_cls1_aprx(:,3));
% hold off;