function [store_concentration]=model(I,k_end,concentration,Gp,RNAP,Ribosome,doubling_time_td,mRNA_halflife,protein_halflife,LX_1,LX_2,LX_3,LT_1,LT_2,LT_3,e_x,e_L,K_IX,K_X,K_IL,K_L,KI1,nI1,W11,WI1,K12,n12,W22,W12,K13,n13,W33,W13,K23,n23,W23)

specific_growth_rate=log(2)/(doubling_time_td);                %1/min
kd_m=log(2)/(mRNA_halflife);                                   %1/min 
kd_p=log(2)/(protein_halflife*60);                             %1/min

KE_1=e_x/LX_1;                                                 %1/min
tau_1=KE_1/K_IX; 
KE_2=e_x/LX_2;                                                 %1/min
tau_2=KE_2/K_IX;
KE_3=e_x/LX_3;                                                 %1/min
tau_3=KE_3/K_IX;    

KE_4=e_L/LT_1;                                                 %1/min
tau_4=KE_4/K_IL;
KE_5=e_L/LT_2;                                                 %1/min                                                                                                       %1/min
tau_5=KE_5/K_IL;  
KE_6=e_L/LT_3;                                                 %1/min                                                                                                       %1/min
tau_6=KE_6/K_IL; 

S=eye(6,6);                 

A=zeros(6,6);                                                  %1/min
    A(1,1)=-(specific_growth_rate+kd_m);
    A(2,2)=-(specific_growth_rate+kd_m);
    A(3,3)=-(specific_growth_rate+kd_m);
    A(4,4)=-(specific_growth_rate+kd_p);
    A(5,5)=-(specific_growth_rate+kd_p);
    A(6,6)=-(specific_growth_rate+kd_p);


step=1.0;                                                      %min

Ahat=expm(A*step);

Shat=A^(-1)*(Ahat-eye(6,6))*S;

store_concentration=zeros(6,k_end+1);

    for k=0:1:k_end

        store_concentration(:,k+1)=concentration;

        fI=(I^nI1)/(KI1^nI1+I^nI1);
        u1=(W11+WI1*fI)/(1+W11+WI1*fI);
        %Specific rate of transcription
        r_x_1=KE_1*RNAP*Gp/(K_X*tau_1+(tau_1+1)*Gp);                                    %nmol/gDWmin
        T_X1=r_x_1*u1;                                                                  %nmol/gDW/min
        %%protein1
        %Specific rate of translation
        r_x_4=KE_4*Ribosome*concentration(1)/(K_L*tau_4+(tau_4+1)*concentration(1));    %nmol/gDWmin
        T_L1=r_x_4;                                                                     %nmol/gDW/min
        %%mRNA2
        %control term
        f12=concentration(4)^n12/(K12^n12+concentration(4)^n12);
        u2=(W22+W12*f12)/(1+W22+W12*f12);                                         
        %Specific rate of transcription
        r_x_2=KE_2*RNAP*Gp/(K_X*tau_2+(tau_2+1)*Gp);                                    %nmol/gDWmin
        T_X2=r_x_2*u2;                                                                  %nmol/gDW/min
        %%protein2
        %Specific rate of translation
        r_x_5=KE_5*Ribosome*concentration(2)/(K_L*tau_5+(tau_5+1)*concentration(2));    %nmol/gDWmin
        T_L2=r_x_5;
        %%mRNA3
        %control term
        f13=concentration(4)^n13/(K13^n13+concentration(4)^n13);
        f23=concentration(5)^n23/(K23^n23+concentration(5)^n23);
        u3=(W33+W13*f13)/(1+W33+W13*f13+W23*f23);                                       %control term
        %Specific rate of transcription
        r_x_3=KE_3*RNAP*Gp/(K_X*tau_3+(tau_3+1)*Gp);                                    %nmol/gDWmin
        T_X3=r_x_3*u3;                                                                  %nmol/gDW/min
        %%protein3
        %Specific rate of translation
        r_x_6=KE_6*Ribosome*concentration(3)/(K_L*tau_6+(tau_6+1)*concentration(3));    %nmol/gDWmiin
        T_L3=r_x_6;                                                                     %nmol/gDW/min

        r=[ T_X1;
            T_X2;
            T_X3;
            T_L1;
            T_L2;
            T_L3];

        new_concentration=Ahat*concentration+Shat*r;
        concentration=new_concentration;
    end

end
