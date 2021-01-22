clear all;  clc;
NN= 420;
MM = 4;
T = 100;
iter = 5000;
betaE = 1;
AN_ratio = 0;
Qtt = 0:5:50;
Qtotal = 10.^(Qtt/10);%10dB
Q0max = AN_ratio*Qtotal;
Qmax = (1-AN_ratio)*Qtotal;
Pmax = 10^0.1;%1dB
sumrate = zeros(length(Qtotal),length(MM));
RS_sum_cls1 = zeros(length(Qtotal),length(MM));
RS_1_cls1 = zeros(length(Qtotal),length(MM));
RS_2_cls1 = zeros(length(Qtotal),length(MM));
RS_3_cls1 = zeros(length(Qtotal),length(MM));
RS_sum_cls1_aprx = zeros(length(Qtotal),length(MM));
RS_1_cls1_aprx = zeros(length(Qtotal),length(MM));
RS_2_cls1_aprx = zeros(length(Qtotal),length(MM));

sumrate_mrt = zeros(length(Qtotal),length(MM));
RS_sum_cls1_mrt = zeros(length(Qtotal),length(MM));
RS_1_cls1_mrt = zeros(length(Qtotal),length(MM));
RS_2_cls1_mrt = zeros(length(Qtotal),length(MM));
RS_sum_cls1_aprx_mrt = zeros(length(Qtotal),length(MM));
RS_1_cls1_aprx_mrt = zeros(length(Qtotal),length(MM));
RS_2_cls1_aprx_mrt = zeros(length(Qtotal),length(MM));


for j=1:length(Q0max)
for i=1:length(MM)
    fprintf('%d of %d\n',i,length(MM));
    PdB = zeros(MM(i),3);
    QdB = zeros(MM(i),3);
    Q0dB = zeros(MM(i),3);
    beta = zeros(MM(i),3);
    for m=1:MM(i)
        PdB(m,1) = Pmax;
        PdB(m,2) = Pmax;
        PdB(m,3) = Pmax;
        QdB(m,1) = Qmax(j)/MM(i)/3;
        QdB(m,2) = Qmax(j)/MM(i)/3;
        QdB(m,3) = Qmax(j)/MM(i)/3;
        Q0dB(m) = Q0max(j)/MM(i);
        beta(m,1) = 1;
        beta(m,2) = 1;
        beta(m,3) = 1;
    end
    tau = MM(i);
    [RSsum,RSsum_aprx,RS1,RS2,RS3,RS1_aprx,RS2_aprx,RS3_aprx,R1,...
        R11,R12,R13,R1_aprx,R11_aprx,R12_aprx,R13_aprx,...
        R2,R2_aprx,R23,R23_aprx,R3,R3_aprx,RE1,RE2,RE3,RE1_aprx,RE2_aprx,RE3_aprx...
        Theta11,I111,I112,I113,I113_sim,I114,I114_sim...
        Theta12,I121,I122,I123,I123_sim,I124,I124_sim...
        Theta13,I131,I132,I133,I133_sim,I134,I134_sim...
        Theta2,Theta23,I21,I22,I23,I231,I232,I233...
        Theta3,I31,I32,I33,H1,H2,H3,Hh,W,z,...
        Q,P,Q0,rho,tmp1,tmp2,tmp3,tmp4,tmp5]...
= Sim_VarN_ZF_3UE(QdB,PdB,Q0dB,beta,betaE,tau,T,MM(i),NN,iter);

[RSsum_mrt,RSsum_aprx_mrt,RS1_mrt,RS2_mrt,RS3_mrt,RS1_aprx_mrt,RS2_aprx_mrt,RS3_aprx_mrt,R1_mrt,...
        R11_mrt,R12_mrt,R13_mrt,R1_aprx_mrt,R11_aprx_mrt,R12_aprx_mrt,R13_aprx_mrt,...
        R2_mrt,R2_aprx_mrt,R23_mrt,R23_aprx_mrt,R3_mrt,R3_aprx_mrt,RE1_mrt,RE2_mrt,RE3_mrt,RE1_aprx_mrt,RE2_aprx_mrt,RE3_aprx_mrt...
        Theta11_mrt,I111_mrt,I112_mrt,I113_mrt,I113_sim_mrt,I114_mrt,I114_sim_mrt...
        Theta12_mrt,I121_mrt,I122_mrt,I123_mrt,I123_sim_mrt,I124_mrt,I124_sim_mrt...
        Theta13_mrt,I131_mrt,I132_mrt,I133_mrt,I133_sim_mrt,I134_mrt,I134_sim_mrt...
        Theta2_mrt,Theta23_mrt,I21_mrt,I22_mrt,I23_mrt,I231_mrt,I232_mrt,I233_mrt...
        Theta3_mrt,I31_mrt,I32_mrt,I33_mrt,H1_mrt,H2_mrt,H3_mrt,Hh_mrt,W_mrt,z_mrt,...
        Q,P,Q0,rho,tmp1,tmp2,tmp3,tmp4,tmp5]...
= Sim_VarN_MRT_3UE(QdB,PdB,Q0dB,beta,betaE,tau,T,MM(i),NN,iter);

sumrate(j,i) = sum(RSsum);
RS_sum_cls1(j,i) = RSsum(:,1);
RS_1_cls1(j,i) = RS1(:,1);
RS_2_cls1(j,i) = RS2(:,1);
RS_3_cls1(j,i) = RS3(:,1);
RS_sum_cls1_aprx(j,i) = RSsum_aprx(:,1);
RS_1_cls1_aprx(j,i) = RS1_aprx(:,1);
RS_2_cls1_aprx(j,i) = RS2_aprx(:,1);

sumrate_mrt(j,i) = sum(RSsum_mrt);
RS_sum_cls1_mrt(j,i) = RSsum_mrt(:,1);
RS_1_cls1_mrt(j,i) = RS1_mrt(:,1);
RS_2_cls1_mrt(j,i) = RS2_mrt(:,1);
RS_sum_cls1_aprx_mrt(j,i) = RSsum_aprx_mrt(:,1);
RS_1_cls1_aprx_mrt(j,i) = RS1_aprx_mrt(:,1);
RS_2_cls1_aprx_mrt(j,i) = RS2_aprx_mrt(:,1);
end



end
% save('MRTvsZF',...
%     'NN','MM','T','betaE',...
%     'AN_ratio','Qtotal','Pmax',...
%     'RS_sum_cls1','RS_1_cls1','RS_2_cls1',...
%     'RS_sum_cls1_aprx','RS_1_cls1_aprx','RS_2_cls1_aprx');
fprintf('DONE!');
figure
% plot(Qtt,RS_sum_cls1(:,1),'Displayname','M=3');
hold on;
% plot(Qtt,RS_sum_cls1(:,2),'Displayname','M=6');
% plot(Qtt,RS_sum_cls1(:,2),'Displayname','M=9');
scatter(Qtt,sumrate);
% scatter(Qtt,RS_sum_cls1_aprx(:,2));
scatter(Qtt,sumrate_mrt);

% plot(Qtt,RS_sum_cls1_mrt(:,1),'Displayname','M=3 mrt');
% % plot(Qtt,RS_sum_cls1_mrt(:,2),'Displayname','M=6 mrt');
% plot(Qtt,RS_sum_cls1_mrt(:,2),'Displayname','M=9 mrt');
% scatter(Qtt,RS_sum_cls1_aprx_mrt(:,1));
% % scatter(Qtt,RS_sum_cls1_aprx_mrt(:,2));
% scatter(Qtt,RS_sum_cls1_aprx_mrt(:,2));
% % hold off;