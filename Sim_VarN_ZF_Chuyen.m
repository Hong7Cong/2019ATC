% simulation & approximation
% calculate RS1, RS2
function[RSsum,RSsum_aprx,RS1,RS2,RS1_aprx,RS2_aprx,R1,...
        R11,R12,R1_aprx,R11_aprx,R12_aprx,...
        R2,R2_aprx,RE1,RE2,RE1_aprx,RE2_aprx,...
        Theta11,I111,I112,I113,I113_sim,I114,I114_sim...
        Theta12,I121,I122,I123,I123_sim,I124,I124_sim...
        H1,H2,Hh,W,z,...
        Q,P,Q0,rho,tmp1,tmp2,tmp3,tmp4,tmp5]...
        = Sim_VarN_ZF_Chuyen(QdB,PdB,Q0dB,beta,betaE,tau,T,MM,NN,iter)

Q = QdB;%power portion for data for each EU
P = PdB;%power of each user in uplink training
Q0 = Q0dB;%power portion for artifical noise for each cluster

rho = zeros(MM,2);%rho_{m,k} as in ref.19-below (4)

Theta11 = zeros(length(NN),MM);
Theta12 = zeros(length(NN),MM);
Theta2 = zeros(length(NN),MM);
ThetaE1 = zeros(length(NN),MM);
ThetaE2 = zeros(length(NN),MM);

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

I21 = zeros(length(NN),MM);
I22 = zeros(length(NN),MM);
I23 = zeros(length(NN),MM);
I22_sim = zeros(length(NN),MM);
I23_sim = zeros(length(NN),MM);

IE11 = zeros(length(NN),MM);
IE12 = zeros(length(NN),MM);
IE13 = zeros(length(NN),MM);

IE21 = zeros(length(NN),MM);
IE22 = zeros(length(NN),MM);
IE23 = zeros(length(NN),MM);


tmp1 = zeros(length(NN),MM);
tmp2 = zeros(length(NN),MM);
tmp3 = zeros(length(NN),MM);
tmp4 = zeros(length(NN),MM);
tmp5 = zeros(length(NN),MM);
tmp6 = zeros(length(NN),MM);
tmp7 = zeros(length(NN),MM);
tmp8 = zeros(length(NN),MM);
tmp9 = zeros(length(NN),MM);

R11 = zeros(length(NN),MM);
R12 = zeros(length(NN),MM);
R1 = zeros(length(NN),MM);
R1_aprx = zeros(length(NN),MM);
R11_aprx = zeros(length(NN),MM);
R12_aprx = zeros(length(NN),MM);
R2 = zeros(length(NN),MM);
R2_aprx = zeros(length(NN),MM);
RE1 = zeros(length(NN),MM);
RE2 = zeros(length(NN),MM);
RE1_aprx = zeros(length(NN),MM);
RE2_aprx = zeros(length(NN),MM);
RS1 = zeros(length(NN),MM);
RS2 = zeros(length(NN),MM);
RS1_aprx = zeros(length(NN),MM);
RS2_aprx = zeros(length(NN),MM);
RSsum = zeros(length(NN),MM);
RSsum_aprx = zeros(length(NN),MM);

InterI1=zeros(length(NN),MM);% Inter-cluster interference (17) for EU1
InterI2=zeros(length(NN),MM);% Inter-cluster interference (17) for EU2

for n=1:length(NN)
    
    %Initialization-------------------------------
    Hh = zeros(NN(n),MM);%channel
    W = zeros(NN(n),MM);%noise
    Y = zeros(NN(n),MM);%received signal at BS
    z = zeros(NN(n),MM);
    
    %For each user in Cluster, calculate rho_{m,k} as in ref.19-below (4)
    for m=1:MM
        rho(m,1) = P(m,1)*beta(m,1)*tau/(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2))*tau);
        rho(m,2) = P(m,2)*beta(m,2)*tau/(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2))*tau);
    end

    %For average------------------------------------
    
    for i=1:iter
        InterI1=InterI1*0;% Reset inter-cluster interference for next loop
        InterI2=InterI2*0;% Reset inter-cluster interference for next loop
        
        H1= sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%complex Gaussian uplink channel from EU1 to BS
        H2= sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%complex Gaussian uplink channel from EU2 to BS
        g = sqrt(0.5)*(randn(NN(n),1) +1j*randn(NN(n),1));%complex Gaussian channel from BS to Eaves.
        N = sqrt(0.5)*(randn(NN(n),MM) +1j*randn(NN(n),MM));%Gaussian noise at BS
        
        %Perfect pilot transmission, channel estimation, and AN caculation
        for m=1:MM %for MM clusters
            Y(:,m) = sqrt(P(m,1)*beta(m,1)*tau)*H1(:,m)...
                + sqrt(P(m,2)*beta(m,2)*tau)*H2(:,m)...
                + N(:,m);%received signal at BS from m-th cluster using (2)
            Hh(:,m) = sqrt((P(m,1)*beta(m,1)+P(m,2)*beta(m,2))*tau)*Y(:,m)...
            /(1+(P(m,1)*beta(m,1)+P(m,2)*beta(m,2))*tau);%MMSE channel estimate for the m-cluster using (3)
           
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
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
                    
                    InterI2(n,m) = InterI2(n,m) + beta(m,2)*...%term (17) for EU2
                        ((Q(j,1)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
                end
            end
            for j=m+1:MM
                InterI1(n,m) = InterI1(n,m) + beta(m,1)*...%term (17) for EU1
                        ((Q(j,1)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
                InterI2(n,m) = InterI2(n,m) + beta(m,2)*...%term (17) for EU2
                        ((Q(j,1)*(abs(ctranspose(H1(:,m))*W(:,j))^2)...
                        +Q(j,2)*(abs(ctranspose(H2(:,m))*W(:,j))^2))...
                        +Q0(j)*(abs(ctranspose(H1(:,m))*z(:,j))^2));
            end
        end        
        
        for m=1:MM
            tmp1(n,m) = tmp1(n,m) + ctranspose(H1(:,m))*W(:,m);%{hm*wm} in (15) for EU1%
            tmp2(n,m) = tmp2(n,m) + abs(ctranspose(H1(:,m))*W(:,m))^2;%{|hm*wm|^2} in (15) for EU1%
            
            tmp3(n,m) = tmp3(n,m) + ctranspose(H2(:,m))*W(:,m);%{hm*wm} in (15) for EU2%
            tmp4(n,m) = tmp4(n,m) + abs(ctranspose(H2(:,m))*W(:,m))^2;%{|hm*wm|^2} in (15) for EU2%
            
            tmp5(n,m) = tmp5(n,m) + abs(ctranspose(g)*W(:,m))^2;%{|g*wm|^2} for Eve%
            
            tmp6(n,m) = tmp6(n,m) + beta(m,1)*... %{hm*zm} in (16) for EU1%
            (Q0(m)*(abs(ctranspose(H1(:,m))*z(:,m))^2));
            tmp7(n,m) = tmp7(n,m) + beta(m,2)*... %{hm*zm} in (16) for EU2%
            (Q0(m)*(abs(ctranspose(H2(:,m))*z(:,m))^2));
        
            tmp8(n,m) = tmp8(n,m) + InterI1(n,m); %term (17) for EU1 Average
            tmp9(n,m) = tmp9(n,m) + InterI2(n,m); %term (17) for EU2 Average
                               
        end
    end
    
    %Average now-----------------------------------------------
    for m=1:MM
    tmp1(n,m) = abs(tmp1(n,m)/iter).^2 ;% |E{hm*wm}|^2 %
    tmp2(n,m) = tmp2(n,m)/iter ;%E{|hm*wm|^2} %
    tmp3(n,m) = abs(tmp3(n,m)/iter).^2 ;% |E{hm*wm}|^2 %
    tmp4(n,m) = tmp4(n,m)/iter ;%E{|hm*wm|^2} %
    tmp5(n,m) = tmp5(n,m)/iter;%E{|g*wm|^2} %
    tmp6(n,m) = tmp6(n,m)/iter;% |E{hm*zm}|^2 %
    tmp7(n,m) = tmp7(n,m)/iter;% |E{hm*zm}|^2 %
    tmp8(n,m) = tmp8(n,m)/iter;% E{term (17)}
    tmp9(n,m) = tmp9(n,m)/iter;% E{term (17)}
    
    %For EU1 to decode itself signal-----------------------------
    Theta11(n,m) = Q(m,1)*beta(m,1)*tmp1(n,m);%desired EU1 signal (13), also second term of (15)
    I111(n,m) = Q(m,1)*beta(m,1)*(tmp2(n,m)-tmp1(n,m));%imperfect channel estimation in (15) for EU1    
    I112(n,m) = Q(m,2)*beta(m,1)*tmp2(n,m);%first term of (16), intra-cluster interference for EU1 
    I113(n,m) = Q0(m)*beta(m,1)*(1-rho(m,1));%second term of (16), AN leakage
    I113_sim(n,m) = tmp6(n,m);%second term of (16), AN leakage, by simulation generating z 
    
    %For EU2 to decode EU1 signal--------------------
    Theta12(n,m) = Q(m,1)*beta(m,2)*tmp3(n,m);%desired EU1 signal (13), also second term of (15)    
    I121(n,m) = Q(m,1)*beta(m,2)*(tmp4(n,m)-tmp3(n,m));%imperfect channel estimation in (15) for EU2
    I122(n,m) = Q(m,2)*beta(m,2)*tmp4(n,m);%first term of (16), intra-cluster interference for EU2
    I123(n,m) = Q0(m)*beta(m,2)*(1-rho(m,2));%second term of (16)
    I123_sim(n,m) = tmp7(n,m);%second term of (16), AN leakage, by simulation generating z
    
    %For EU2 to decode itself signal------------------
    Theta2(n,m) = Q(m,2)*beta(m,2)*tmp3(n,m);%desired EU2 signal (13) 
    I21(n,m) = Q(m,2)*beta(m,2)*(tmp4(n,m)-tmp3(n,m));%imperfect channel estimation in (15) for EU2
    I22(n,m) = Q0(m)*beta(m,2)*(1-rho(m,2));%second term of (16), AN leakage
    I22_sim(n,m) = tmp7(n,m); %second term of (16), AN leakage , by simulation generating z
    
    %For Eve------------------------------------------------------
    ThetaE1(n,m) = Q(m,1)*betaE*tmp5(n,m);%desired Eve. signal (E1)
    IE11(n,m) = Q(m,2)*betaE*tmp5(n,m);%first term of E2
    IE12(n,m) = Q0(m)*betaE;%second term of E2
    
    ThetaE2(n,m) = Q(m,2)*betaE*tmp5(n,m);
    IE21(n,m) = Q(m,1)*betaE*tmp5(n,m);
    IE22(n,m) = Q0(m)*betaE;%second term of E2
    
    %For inter-cluster interference in (17) and E3------------------------
    if(m>1)    
        for j=1:m-1
            I114(n,m) = I114(n,m) + beta(m,1)*((1-rho(m,1))*(Q(j,1)+Q(j,2))+Q0(j));
            I124(n,m) = I124(n,m) + beta(m,2)*((1-rho(m,2))*(Q(j,1)+Q(j,2))+Q0(j));
            IE13(n,m) = IE13(n,m) + betaE*(Q(j,1)+Q(j,2)+Q0(j)); 
            IE23(n,m) = IE13(n,m);
        end
    end
    for j=m+1:MM
        I114(n,m) = I114(n,m) + beta(m,1)*((1-rho(m,1))*(Q(j,1)+Q(j,2))+Q0(j));
        I124(n,m) = I124(n,m) + beta(m,2)*((1-rho(m,2))*(Q(j,1)+Q(j,2))+Q0(j));
        IE13(n,m) = IE13(n,m) + betaE*(Q(j,1)+Q(j,2)+Q0(j));
        IE23(n,m) = IE13(n,m);
    end
    I114_sim(n,m) = tmp8(n,m);
    I124_sim(n,m) = tmp9(n,m);
    I23(n,m) = I124(n,m);%used with theta2
    I23_sim(n,m) = I124_sim(n,m);%used with theta2
    
    %rate at EU1---------------------------------------
    R11(n,m) = (1-tau/T)*log2(1+Theta11(n,m)/...
        (I111(n,m)+I112(n,m)+I113_sim(n,m)+I114_sim(n,m)+1));%due to EU1 decode itself
    R12(n,m) = (1-tau/T)*log2(1+Theta12(n,m)/...
        (I121(n,m)+I122(n,m)+I123_sim(n,m)+I124_sim(n,m)+1));%due to EU2 decode EU1
    R1(n,m) = min(R11(n,m),R12(n,m));        
    
    %rate at EU1 approximation-------------------------
    R11_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*beta(m,1)*rho(m,1)*(NN(n)-MM+1)/...
        (I111(n,m)+I112(n,m)+I113(n,m)+I114(n,m)+1));%due to EU1 decode itself
    R12_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*beta(m,2)*rho(m,2)*(NN(n)-MM+1)/...
        (I121(n,m)+I122(n,m)+I123(n,m)+I124(n,m)+1));%due to EU2 decode EU1         
    R1_aprx(n,m) = min(R12_aprx(n,m),R11_aprx(n,m));
    
    %rate at EU2---------------------------------------
    R2(n,m) = (1-tau/T)*log2(1+Theta2(n,m)/...
        (I21(n,m)+I22_sim(n,m)+I23_sim(n,m)+1));        
    R2_aprx(n,m) = (1-tau/T)*log2(1+Q(m,2)*beta(m,2)*rho(m,2)*(NN(n)-MM+1)/...
        (I21(n,m)+I22(n,m)+I23(n,m)+1));%need to recalculate 
    
    %rate at Eve---------------------------------------
    RE1(n,m) = (1-tau/T)*log2(1+ThetaE1(n,m))/...
        (IE11(n,m) + IE12(n,m) + IE13(n,m));
    RE2(n,m) = (1-tau/T)*log2(1+ThetaE2(n,m))/...
        (IE21(n,m) + IE22(n,m) + IE23(n,m));
    RE1_aprx(n,m) = (1-tau/T)*log2(1+Q(m,1)*betaE)/...
        (IE11(n,m) + IE12(n,m) + IE13(n,m));
    RE2_aprx(n,m) = (1-tau/T)*log2(1+Q(m,2)*betaE)/...
        (IE21(n,m) + IE22(n,m) + IE23(n,m));
    
    %secrecy rate ----------------------
    RS1(n,m) = max(0, R1(n,m)-RE1(n,m));
    RS2(n,m) = max(0, R2(n,m)-RE2(n,m));        
    RS1_aprx(n,m) = max(0, R1_aprx(n,m)-RE1_aprx(n,m));
    RS2_aprx(n,m) = max(0, R2_aprx(n,m)-RE2_aprx(n,m));
    
    %sumrate
    RSsum(n,m) =  RS1(n,m)+ RS2(n,m);
    RSsum_aprx(n,m) =  RS1_aprx(n,m)+ RS2_aprx(n,m);
    end    
end 

end