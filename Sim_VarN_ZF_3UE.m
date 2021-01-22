% simulation & approximation
% calculate RS1, RS2
function[RSsum,RSsum_aprx,RS1,RS2,RS3,RS1_aprx,RS2_aprx,RS3_aprx,R1,...
        R11,R12,R13,R1_aprx,R11_aprx,R12_aprx,R13_aprx,...
        R2,R2_aprx,R23,R23_aprx,R3,R3_aprx,RE1,RE2,RE3,RE1_aprx,RE2_aprx,RE3_aprx...
        Theta11,I111,I112,I113,I113_sim,I114,I114_sim...
        Theta12,I121,I122,I123,I123_sim,I124,I124_sim...
        Theta13,I131,I132,I133,I133_sim,I134,I134_sim...
        Theta2,Theta23,I21,I22,I23,I231,I232,I233...
        Theta3,I31,I32,I33,H1,H2,H3,Hh,W,z,...
        Q,P,Q0,rho,tmp1,tmp2,tmp3,tmp4,tmp5]...
        = Sim_VarN_ZF_3UE(QdB,PdB,Q0dB,beta,betaE,tau,T,MM,NN,iter)

Q = QdB;%power portion for data for each EU
P = PdB;%power of each user in uplink training
Q0 = Q0dB;%power portion for artifical noise for each cluster

rho = zeros(MM,2);%rho_{m,k} as in ref.19-below (4)

Theta11 = zeros(length(NN),MM);
Theta12 = zeros(length(NN),MM);
Theta13 = zeros(length(NN),MM);
Theta2 = zeros(length(NN),MM);
Theta23 = zeros(length(NN),MM);
Theta3 = zeros(length(NN),MM);
ThetaE1 = zeros(length(NN),MM);
ThetaE2 = zeros(length(NN),MM);
ThetaE3 = zeros(length(NN),MM);

I111 = zeros(length(NN),MM);
I112 = zeros(length(NN),MM);
I113 = zeros(length(NN),MM);
I113_sim = zeros(length(NN),MM);
I114 = zeros(length(NN),MM);
I114_sim = zeros(length(NN),MM);

I121 = zeros(length(NN),MM);
I122 = zeros(length(NN),MM);
I123 = zeros(length(NN),MM);
I123_sim = zeros(length(NN),MM);
I124 = zeros(length(NN),MM);
I124_sim = zeros(length(NN),MM);

I131 = zeros(length(NN),MM);
I132 = zeros(length(NN),MM);
I133 = zeros(length(NN),MM);
I133_sim = zeros(length(NN),MM);
I134 = zeros(length(NN),MM);
I134_sim = zeros(length(NN),MM);

I21 = zeros(length(NN),MM);
I22 = zeros(length(NN),MM);
I23 = zeros(length(NN),MM);
I22_sim = zeros(length(NN),MM);
I23_sim = zeros(length(NN),MM);

I231 = zeros(length(NN),MM);
I232 = zeros(length(NN),MM);
I233 = zeros(length(NN),MM);
I232_sim = zeros(length(NN),MM);
I233_sim = zeros(length(NN),MM);

I31 = zeros(length(NN),MM);
I32 = zeros(length(NN),MM);
I33 = zeros(length(NN),MM);
I32_sim = zeros(length(NN),MM);
I33_sim = zeros(length(NN),MM);

IE11 = zeros(length(NN),MM);
IE12 = zeros(length(NN),MM);
IE13 = zeros(length(NN),MM);

IE21 = zeros(length(NN),MM);
IE22 = zeros(length(NN),MM);
IE23 = zeros(length(NN),MM);

IE31 = zeros(length(NN),MM);
IE32 = zeros(length(NN),MM);
IE33 = zeros(length(NN),MM);

tmp1 = zeros(length(NN),MM);
tmp2 = zeros(length(NN),MM);
tmp3 = zeros(length(NN),MM);
tmp4 = zeros(length(NN),MM);
tmp5 = zeros(length(NN),MM);
tmp6 = zeros(length(NN),MM);
tmp7 = zeros(length(NN),MM);
tmp8 = zeros(length(NN),MM);
tmp9 = zeros(length(NN),MM);
tmp10 = zeros(length(NN),MM);%ADD
tmp11 = zeros(length(NN),MM);%ADD
tmp12 = zeros(length(NN),MM);%ADD
tmp13 = zeros(length(NN),MM);%ADD

R11 = zeros(length(NN),MM);
R12 = zeros(length(NN),MM);
R13 = zeros(length(NN),MM);%ADD

R1 = zeros(length(NN),MM);
R1_aprx = zeros(length(NN),MM);
R11_aprx = zeros(length(NN),MM);
R12_aprx = zeros(length(NN),MM);
R13_aprx = zeros(length(NN),MM);
R2 = zeros(length(NN),MM);
R22 = zeros(length(NN),MM);
R22_aprx = zeros(length(NN),MM);
R2_aprx = zeros(length(NN),MM);
R23 = zeros(length(NN),MM);
R23_aprx = zeros(length(NN),MM);
R3 = zeros(length(NN),MM);
R3_aprx = zeros(length(NN),MM);

RE1 = zeros(length(NN),MM);
RE2 = zeros(length(NN),MM);
RE3 = zeros(length(NN),MM);
RE1_aprx = zeros(length(NN),MM);
RE2_aprx = zeros(length(NN),MM);
RE3_aprx = zeros(length(NN),MM);

RS1 = zeros(length(NN),MM);
RS2 = zeros(length(NN),MM);
RS3 = zeros(length(NN),MM);
RS1_aprx = zeros(length(NN),MM);
RS2_aprx = zeros(length(NN),MM);
RS3_aprx = zeros(length(NN),MM);
RSsum = zeros(length(NN),MM);
RSsum_aprx = zeros(length(NN),MM);

InterI1=zeros(length(NN),MM);% Inter-cluster interference (17) for EU1
InterI2=zeros(length(NN),MM);% Inter-cluster interference (17) for EU2
InterI3=zeros(length(NN),MM);% Inter-cluster interference (17) for EU3

for n=1:length(NN)
    
    %Initialization-------------------------------
    Hh = zeros(NN(n),MM);%channel
    W = zeros(NN(n),MM);%noise
    Y = zeros(NN(n),MM);%received signal at BS
    z = zeros(NN(n),MM);
    
    %For each user in Cluster, calculate rho_{m,k} as in ref.19-below (4)
    for m=1:MM
        rho(m,1) = P(m,1)*beta(m,1)*tau/(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2)+P(m,3)*beta(m,3))*tau);%ADD beta(m,3)
        rho(m,2) = P(m,2)*beta(m,2)*tau/(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2)+P(m,3)*beta(m,3))*tau);%ADD beta(m,3)
        rho(m,3) = P(m,3)*beta(m,3)*tau/(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2)+P(m,3)*beta(m,3))*tau);%ADD beta(m,3)
    end

    %For average------------------------------------
    
    for i=1:iter
        InterI1=InterI1*0;% Reset inter-cluster interference for next loop
        InterI2=InterI2*0;% Reset inter-cluster interference for next loop
        InterI3=InterI3*0;% Reset inter-cluster interference for next loop
        
        H1= sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%complex Gaussian uplink channel from EU1 to BS
        H2= sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%complex Gaussian uplink channel from EU2 to BS
        H3= sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%complex Gaussian uplink channel from EU2 to BS %ADD
        g = sqrt(0.5)*(randn(NN(n),1) +1j*randn(NN(n),1));%complex Gaussian channel from BS to Eaves.
        N = sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%Gaussian noise at BS
        
        %Perfect pilot transmission, channel estimation, and AN caculation
        for m=1:MM %for MM clusters
            Y(:,m) = sqrt(P(m,1)*beta(m,1)*tau)*H1(:,m)...
                + sqrt(P(m,2)*beta(m,2)*tau)*H2(:,m)...
                + sqrt(P(m,3)*beta(m,3)*tau)*H3(:,m)...%ADD
                + N(:,m);%received signal at BS from m-th cluster using (2)
            Hh(:,m) = sqrt((P(m,1)*beta(m,1)+P(m,2)*beta(m,2)+P(m,3)*beta(m,3))*tau)*Y(:,m)...%ADD P(m,3)*beta(m,3)
            /(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2)+P(m,3)*beta(m,3))*tau);%MMSE channel estimate for the m-cluster using (3)
           
%             % Calculate AN
%             temp = null(ctranspose(Hh(:,m)));
%             z(:,m) = temp(:,1);
        end
        %Calculate AN, which is in the null space of channel estimate
        z= null(ctranspose(Hh));
        
        %Precoding matrix using (5)
        V = Hh*(ctranspose(Hh)*Hh)^(-1);
        for m=1:MM
%             W(:,m) = Hh(:,m)/norm(Hh(:,m));
            W(:,m) = V(:,m)/norm(V(:,m));
        end
        
        %--------Calculate Intercell Interference by simulation (17)--------------
        

        for m=1:MM
            if(m>1)    
                for j=1:m-1
                    InterI1(n,m) = InterI1(n,m) + beta(m,1)*...%term (17) for EU1
                        ((Q(j,1)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H1(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
                    
                    InterI2(n,m) = InterI2(n,m) + beta(m,2)*...%term (17) for EU2
                        ((Q(j,1)*(abs(ctranspose(H2(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H2(:,m))*z(:,j))^2));
                    
                    InterI3(n,m) = InterI3(n,m) + beta(m,3)*...%term (17) for EU3
                        ((Q(j,1)*(abs(ctranspose(H3(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H3(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H3(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H3(:,m))*z(:,j))^2));
                end
            end
            for j=m+1:MM
                InterI1(n,m) = InterI1(n,m) + beta(m,1)*...%term (17) for EU1
                        ((Q(j,1)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H1(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
                InterI2(n,m) = InterI2(n,m) + beta(m,2)*...%term (17) for EU2
                        ((Q(j,1)*(abs(ctranspose(H2(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H2(:,m))*z(:,j))^2));
                InterI3(n,m) = InterI3(n,m) + beta(m,2)*...%term (17) for EU3
                        ((Q(j,1)*(abs(ctranspose(H3(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H3(:,m))*W(:,j))^2)...
                        +Q(j,3)*(abs(ctranspose(H3(:,m))*W(:,j))^2))...%ADD
                        +Q0(j)*(abs(ctranspose(H3(:,m))*z(:,j))^2));
            end
        end        
        
        for m=1:MM
            tmp1(n,m) = tmp1(n,m) + ctranspose(H1(:,m))*W(:,m);%{hm*wm} in (15) for EU1%
            tmp2(n,m) = tmp2(n,m) + abs(ctranspose(H1(:,m))*W(:,m))^2;%{|hm*wm|^2} in (15) for EU1%            
            tmp3(n,m) = tmp3(n,m) + ctranspose(H2(:,m))*W(:,m);%{hm*wm} in (15) for EU2%
            tmp4(n,m) = tmp4(n,m) + abs(ctranspose(H2(:,m))*W(:,m))^2;%{|hm*wm|^2} in (15) for EU2%
            tmp10(n,m) = tmp10(n,m) + ctranspose(H3(:,m))*W(:,m);%{hm*wm} in (15) for EU3% ADD
            tmp11(n,m) = tmp11(n,m) + abs(ctranspose(H3(:,m))*W(:,m))^2;%{|hm*wm|^2} in (15) for EU3% ADD
            
            tmp5(n,m) = tmp5(n,m) + abs(ctranspose(g)*W(:,m))^2;%{|g*wm|^2} for Eve%
            
            tmp6(n,m) = tmp6(n,m) + beta(m,1)*... %{hm*zm} in (16) for EU1%
            (Q0(m)*(abs(ctranspose(H1(:,m))*z(:,m))^2));
            tmp7(n,m) = tmp7(n,m) + beta(m,2)*... %{hm*zm} in (16) for EU2%
            (Q0(m)*(abs(ctranspose(H2(:,m))*z(:,m))^2));
            tmp12(n,m) = tmp12(n,m) + beta(m,3)*... %{hm*zm} in (16) for EU3%
            (Q0(m)*(abs(ctranspose(H3(:,m))*z(:,m))^2));
        
            tmp8(n,m) = tmp8(n,m) + InterI1(n,m); %term (17) for EU1 Average
            tmp9(n,m) = tmp9(n,m) + InterI2(n,m); %term (17) for EU2 Average
            tmp13(n,m)=tmp13(n,m) + InterI3(n,m); %term (17) for EU3 Average                
        end
    end
    
    %Average now-----------------------------------------------
    for m=1:MM
    tmp1(n,m) = abs(tmp1(n,m)/iter).^2 ;% |E{hm*wm}|^2 %          EU1
    tmp2(n,m) = tmp2(n,m)/iter ;%E{|hm*wm|^2} %                   EU1
    tmp3(n,m) = abs(tmp3(n,m)/iter).^2 ;% |E{hm*wm}|^2 %          EU2
    tmp4(n,m) = tmp4(n,m)/iter ;%E{|hm*wm|^2} %                   EU2
    tmp10(n,m) = abs(tmp10(n,m)/iter).^2 ;% |E{hm*wm}|^2 %        EU3
    tmp11(n,m) = tmp11(n,m)/iter ;%E{|hm*wm|^2} %                 EU3
    
    tmp5(n,m) = tmp5(n,m)/iter;%E{|g*wm|^2} %
    tmp6(n,m) = tmp6(n,m)/iter;% |E{hm*zm}|^2 %
    tmp7(n,m) = tmp7(n,m)/iter;% |E{hm*zm}|^2 %
    tmp12(n,m) = tmp12(n,m)/iter;% |E{hm*zm}|^2 %
    tmp8(n,m) = tmp8(n,m)/iter;% E{term (17)} EU1
    tmp9(n,m) = tmp9(n,m)/iter;% E{term (17)} EU2
    tmp13(n,m) = tmp13(n,m)/iter;% E{term (17)} EU3
    
    %For EU1 to decode itself signal-----------------------------
    Theta11(n,m) = Q(m,1)*beta(m,1)*tmp1(n,m);%desired EU1 signal (13), also second term of (15)
    I111(n,m) = Q(m,1)*beta(m,1)*(tmp2(n,m)-tmp1(n,m));%imperfect channel estimation in (15) for EU1    
    I112(n,m) = (Q(m,2)+Q(m,3))*beta(m,1)*tmp2(n,m);%first term of (16), intra-cluster interference for EU1 from UE2 UE3
    I113(n,m) = 0;%second term of (16), AN leakage
    I113_sim(n,m) = 0;%second term of (16), AN leakage, by simulation generating z 
    %For EU2 to decode EU1 signal--------------------
    Theta12(n,m) = Q(m,1)*beta(m,2)*tmp3(n,m);%desired EU1 signal (13), also second term of (15)    
    I121(n,m) = Q(m,1)*beta(m,2)*(tmp4(n,m)-tmp3(n,m));%imperfect channel estimation in (15) for EU2
    I122(n,m) = (Q(m,2)+Q(m,3))*beta(m,2)*tmp4(n,m);%first term of (16), intra-cluster interference for EU2
    I123(n,m) = 0;%second term of (16)
    I123_sim(n,m) = 0;%second term of (16), AN leakage, by simulation generating z
    %For EU3 to decode EU1 signal--------------------
    Theta13(n,m) = Q(m,1)*beta(m,3)*tmp10(n,m);%desired EU1 signal (13), also second term of (15)    
    I131(n,m) = Q(m,1)*beta(m,3)*(tmp11(n,m)-tmp10(n,m));%imperfect channel estimation in (15) for EU3
    I132(n,m) = (Q(m,2)+Q(m,3))*beta(m,3)*tmp11(n,m);%first term of (16), intra-cluster interference 
    I133(n,m) = 0;%second term of (16)
    I133_sim(n,m) = 0;%second term of (16), AN leakage, by simulation generating z
    
    %For EU2 to decode itself signal------------------
    Theta2(n,m) = Q(m,2)*beta(m,2)*tmp3(n,m);%desired EU2 signal (13) 
    I21(n,m) = Q(m,2)*beta(m,2)*(tmp4(n,m)-tmp3(n,m));%imperfect channel estimation in (15) for EU2
    I22(n,m) = Q(m,3)*beta(m,3)*tmp4(n,m);%first term of (16), intra-cluster interference from UE3 only because SIC extract interference from UE1
    I22_sim(n,m) = 0; %second term of (16), AN leakage , by simulation generating z   
    %For EU3 to decode EU2 signal--------------------
    Theta23(n,m) = Q(m,2)*beta(m,3)*tmp10(n,m);%desired EU2 signal (13) 
    I231(n,m) = Q(m,2)*beta(m,3)*(tmp11(n,m)-tmp10(n,m));%imperfect channel estimation in (15) for EU2
    I232(n,m) = Q(m,3)*beta(m,3)*tmp11(n,m);%first term of (16), intra-cluster interference from UE3 only because SIC extract interference from UE1
    I232_sim(n,m) = 0; %second term of (16), AN leakage , by simulation generating z
    
    
    %For EU3 to decode itself signal------------------
    Theta3(n,m) = Q(m,3)*beta(m,3)*tmp10(n,m);%desired EU2 signal (13) 
    I31(n,m) = Q(m,3)*beta(m,3)*(tmp11(n,m)-tmp10(n,m));%imperfect channel estimation in (15) for EU3
    I32(n,m) = 0;%first term of (16), intra-cluster interference from UE3 after SIC extract UE1 UE2
    I32_sim(n,m) = 0; %second term of (16), AN leakage , by simulation generating z
    
    %For Eve------------------------------------------------------
    ThetaE1(n,m) = Q(m,1)*betaE*tmp5(n,m);%desired Eve. signal (E1)
    IE11(n,m) = (Q(m,2)+Q(m,3))*betaE*tmp5(n,m);%first term of E2
    IE12(n,m) = 0;%second term of E2
    
    ThetaE2(n,m) = Q(m,2)*betaE*tmp5(n,m);
    IE21(n,m) = (Q(m,1)+Q(m,3))*betaE*tmp5(n,m);
    IE22(n,m) = 0;%second term of E2
    
    ThetaE3(n,m) = Q(m,3)*betaE*tmp5(n,m);
    IE31(n,m) = (Q(m,1)+Q(m,2))*betaE*tmp5(n,m);
    IE32(n,m) = 0;
    
    %For inter-cluster interference in (17) and E3------------------------
    if(m>1)    
        for j=1:m-1
            I114(n,m) = I114(n,m) + beta(m,1)*(1-rho(m,1))*(Q(j,1)+Q(j,2)+Q(j,3));%ADD Q(j,3)
            I124(n,m) = I124(n,m) + beta(m,2)*(1-rho(m,2))*(Q(j,1)+Q(j,2)+Q(j,3));%ADD Q(j,3)
            I134(n,m) = I134(n,m) + beta(m,3)*(1-rho(m,3))*(Q(j,1)+Q(j,2)+Q(j,3));%ADD
            IE13(n,m) = IE13(n,m) + betaE*(Q(j,1)+Q(j,2)+Q(j,3)+Q0(j)); %ADD Q(j,3)
            IE23(n,m) = IE13(n,m);
            IE33(n,m) = IE13(n,m);%ADD
        end
    end
    for j=m+1:MM
        I114(n,m) = I114(n,m) + beta(m,1)*((1-rho(m,1))*(Q(j,1)+Q(j,2)+Q(j,3))+Q0(j));%ADD Q(j,3)
        I124(n,m) = I124(n,m) + beta(m,2)*((1-rho(m,2))*(Q(j,1)+Q(j,2)+Q(j,3))+Q0(j));%ADD Q(j,3)
        I134(n,m) = I134(n,m) + beta(m,3)*((1-rho(m,3))*(Q(j,1)+Q(j,2)+Q(j,3))+Q0(j));%ADD
        IE13(n,m) = IE13(n,m) + betaE*(Q(j,1)+Q(j,2)+Q0(j)); %ADD Q(j,3)
        IE23(n,m) = IE13(n,m);
        IE33(n,m) = IE13(n,m);%ADD
    end
    I114_sim(n,m) = tmp8(n,m);
    I124_sim(n,m) = tmp9(n,m);
    I23(n,m) = I124(n,m);%used with theta2
    I23_sim(n,m) = I124_sim(n,m);%used with theta2
    I33(n,m) = I134(n,m);%add
    I233(n,m) = I124(n,m);%add
    %rate at EU1---------------------------------------
    R11(n,m) = (1-tau/T)*log2(1+Theta11(n,m)/...
        (I111(n,m)+I112(n,m)+I113(n,m)+I114(n,m)+1));%due to EU1 decode itself
    R12(n,m) = (1-tau/T)*log2(1+Theta12(n,m)/...
        (I121(n,m)+I122(n,m)+I123(n,m)+I124(n,m)+1));%due to EU2 decode EU1
    R13(n,m) = (1-tau/T)*log2(1+Theta13(n,m)/...
        (I131(n,m)+I132(n,m)+I133(n,m)+I134(n,m)+1));%due to EU3 decode EU1
    R1(n,m) = min([R11(n,m) R12(n,m) R13(n,m)]);        
    
    %rate at EU1 approximation-------------------------
    R11_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*beta(m,1)*rho(m,1)*(NN(n)-MM+1)/...
        (I111(n,m)+I112(n,m)+I113(n,m)+I114(n,m)+1));%due to EU1 decode itself
    R12_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*beta(m,2)*rho(m,2)*(NN(n)-MM+1)/...
        (I121(n,m)+I122(n,m)+I123(n,m)+I124(n,m)+1));%due to EU2 decode EU1      
    R13_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*beta(m,3)*rho(m,3)*(NN(n)-MM+1)/...
        (I131(n,m)+I132(n,m)+I133(n,m)+I134(n,m)+1));%due to EU3 decode EU1    
    R1_aprx(n,m) = min([R12_aprx(n,m) R11_aprx(n,m) R13_aprx(n,m)]);
    
    %rate at EU2---------------------------------------
    R22(n,m) = (1-tau/T)*log2(1+Theta2(n,m)/...%due to EU2 decode itself
        (I21(n,m)+I22(n,m)+I23(n,m)+1));        
    R22_aprx(n,m) = (1-tau/T)*log2(1+Q(m,2)*beta(m,2)*rho(m,2)*(NN(n)-MM+1)/...
        (I21(n,m)+I22(n,m)+I23(n,m)+1));
    R23(n,m) = (1-tau/T)*log2(1+Theta23(n,m)/...%due to EU3 decode UE2
        (I231(n,m)+I232(n,m)+I233(n,m)+1));     
    R23_aprx(n,m) = (1-tau/T)*log2(1+Q(m,3)*beta(m,3)*rho(m,3)*(NN(n)-MM+1)/...
        (I231(n,m)+I232(n,m)+I233(n,m)+1));
    R2(n,m) = min(R22(n,m),R23(n,m));  
    R2_aprx(n,m) = min(R22_aprx(n,m),R23_aprx(n,m));  
    
    %rate at EU3---------------------------------------
    R3(n,m) = (1-tau/T)*log2(1+Theta3(n,m)/...
    (I31(n,m)+I32(n,m)+I33(n,m)+1));        
    R3_aprx(n,m) = (1-tau/T)*log2(1+Q(m,3)*beta(m,3)*rho(m,3)*(NN(n)-MM+1)/...
    (I31(n,m)+I32(n,m)+I33(n,m)+1));%need to recalculate 
    
    %rate at Eve---------------------------------------
    RE1(n,m) = (1-tau/T)*log2(1+ThetaE1(n,m))/...
        (IE11(n,m) + IE12(n,m) + IE13(n,m));
    RE2(n,m) = (1-tau/T)*log2(1+ThetaE2(n,m))/...
        (IE21(n,m) + IE22(n,m) + IE23(n,m));
    RE3(n,m) = (1-tau/T)*log2(1+ThetaE3(n,m))/...
        (IE31(n,m) + IE32(n,m) + IE33(n,m));
    RE1_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*betaE)/...
        (IE11(n,m) + IE12(n,m) + IE13(n,m));
    RE2_aprx(n,m) = (1-tau/T)*log2(1+Q(m,2)*betaE)/...
        (IE21(n,m) + IE22(n,m) + IE23(n,m));
    RE3_aprx(n,m) = (1-tau/T)*log2(1+Q(m,3)*betaE)/...
        (IE31(n,m) + IE32(n,m) + IE33(n,m));
    
    %secrecy rate ----------------------
    RS1(n,m) = max(0, R1(n,m)-RE1(n,m));
    RS2(n,m) = max(0, R2(n,m)-RE2(n,m));       
    RS3(n,m) = max(0, R3(n,m)-RE3(n,m));  
    
    RS1_aprx(n,m) = max(0, R1_aprx(n,m)-RE1_aprx(n,m));
    RS2_aprx(n,m) = max(0, R2_aprx(n,m)-RE2_aprx(n,m));
    RS3_aprx(n,m) = max(0, R3_aprx(n,m)-RE3_aprx(n,m));
    %sumrate
    RSsum(n,m) =  RS1(n,m)+ RS2(n,m)+ RS3(n,m);
    RSsum_aprx(n,m) =  RS1_aprx(n,m)+ RS2_aprx(n,m)+ RS3_aprx(n,m);
    end    
end 

end