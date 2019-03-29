%Samavi Farnush Bint E Naser
%CHEME 7770
%Prelim#1
%Q2
%21 Mar 2019
%--------------------------------------------------------------------------%

clc
clear all

%%part(a)
%--------------------------------------------------------------------------%

%%general parameters
AvN=6.023*10^23;                                           %Avegadro number
cell_weight=2.8*10^-13;                                    %g/cell 
water_fraction=0.70;
conversion_factor=10^3/(cell_weight*(1-water_fraction));   %uM==>nmol/gDW
RNAP=1150*10^6/AvN*conversion_factor;                      %nmol/gDW 
G=200*10^6/(AvN)*conversion_factor;                        %nmol/gDW
Ribosome=45000*10^6/(AvN)*conversion_factor;               %nmol/gDW 
doubling_time_td=40;                                       %min
mRNA_halflife=2.1;                                         %min
protein_halflife=24;                                       %hr 

%read length
LX_1=1200;                                                 %nt
LT_1=1200/3;                                               %aa
LX_2=2400;                                                 %nt
LT_2=2400/3;                                               %aa
LX_3=600;                                                  %nt
LT_3=600/3;                                                %aa

e_x=60*60;                                                 %nt/min 
e_L=16.5*60;                                               %aa/min 
transcription_initiation=60/42;                            %McClure min^-1                                                    
saturation_transcription=0.24;                             %nmol/gDW 
translation_initiation=60/15;                              %BNID:109525 1/min
saturation_translation=454.64;                             %nmol/gDW 

%promoter parameters
%P1
W11=0.0000001;
KI1=0.3;                                                   %mM
nI1=1.5;                                  
WI1=100;     
%P2
W22=0.0000001;
K12=1000.0;                                                %nmol/gDW
n12=1.5;                             
W12=10;    
%P3
W33=0.0000001;
K13=1000;                                                  %nmol/gDW
n13=1.5;
W13=5;
K23=10000.0;                                               %nmol/gDW
n23=10;
W23=50;                          

parameter = [
    G,RNAP,Ribosome,doubling_time_td,mRNA_halflife,protein_halflife,LX_1,LX_2,LX_3,LT_1,LT_2,LT_3,e_x,e_L,transcription_initiation,saturation_transcription,translation_initiation,saturation_translation,KI1,nI1,W11,WI1,K12,n12,W22,W12,K13,n13,W33,W13,K23,n23,W23
    ];

%get to steady state without inducer
initial_concentration=[
                        0.0;                               %mRNA1
                        0.0;                               %mRNA2
                        0.0;                               %mRNA3
                        0.0;                               %protein1
                        0.0;                               %protein2
                        0.0;                               %protein3
      ];
end_sim=350;
inducer=0;
ss_time=[0:end_sim];
[ss_concentration]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));

figure(1)
hold on
plot(ss_time,ss_concentration(4,:),'g-',ss_time,ss_concentration(5,:),'b-',ss_time,ss_concentration(6,:),'k-');
legend('protein 1','protein 2','protein 3')
xlabel("time (min)")
ylabel("Concentration (nmol/gDW)")
title("steady state without inducer")
legend('protein1','protein 2','protein 3')

%run for an additional 60 min after reaching steady state_phase 1
initial_concentration=ss_concentration(:,end);
end_sim=60;
inducer=0;
time1=[0:end_sim];
[concentration1]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));

%run for 300 min with inducer = 10 mM_phase 2
initial_concentration=concentration1(:,end);
end_sim=300;
inducer=10;
time2=[time1(end):time1(end)+end_sim];
[concentration2]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));

time=[time1 time2];
concentration=[concentration1 concentration2];

figure(2)
hold on
plot(time,concentration(4,:),'g-',time,concentration(5,:),'b-',time,concentration(6,:),'k-');
legend('protein 1','protein 2','protein 3')
xlabel("time (min)")
ylabel("Concentration (nmol/gDW)")
xlim([0 360]);
legend('protein1','protein 2','protein 3')
title("protein evolution with time")

%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
%%part(b)
%--------------------------------------------------------------------------%

save_conc=[concentration1 concentration2];
save_param=parameter;
[list_state,list_datapoint]=size(save_conc);

%get modified parameters
%to implement central difference scheme I am choosing to get both up and
%down regulated parameter arrays

list_param=length(parameter);
delta1=0.05;
parameter_array_up=getparam(list_param,parameter,delta1);
delta2=-0.05;
parameter_array_down=getparam(list_param,parameter,delta2);
param_array=[parameter_array_up;parameter_array_down];

%collect data for modified parameters
%create collection array
phase1=[];
phase2=[];
phase3=[];

for index=1:length(param_array)
    parameter=param_array(index,:);
    %get to steady state without inducer (to ensure steady state with
    %modifed parameters run simulaion for a long time
    initial_concentration=[
                        0.0;                               %mRNA1
                        0.0;                               %mRNA2
                        0.0;                               %mRNA3
                        0.0;                               %protein1
                        0.0;                               %protein2
                        0.0;                               %protein3
                        ];
     end_sim=600;
     inducer=0;
     ss_time=[0:end_sim];
     [ss_concentration]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));
     
     %run for an additional 60 min after reaching steady state_phase 1
     initial_concentration=ss_concentration(:,end);
     end_sim=60;
     inducer=0;
     time1=[0:end_sim];
     [concentration1]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));

     %run for 300 min with inducer = 10 mM_phase 2
     initial_concentration=concentration1(:,end);
     end_sim=300;
     inducer=10;
     time2=[time1(end):time1(end)+end_sim];
     [concentration2]=model(inducer,end_sim,initial_concentration,parameter(1),parameter(2),parameter(3),parameter(4),parameter(5),parameter(6),parameter(7),parameter(8),parameter(9),parameter(10),parameter(11),parameter(12),parameter(13),parameter(14),parameter(15),parameter(16),parameter(17),parameter(18),parameter(19),parameter(20),parameter(21),parameter(22),parameter(23),parameter(24),parameter(25),parameter(26),parameter(27),parameter(28),parameter(29),parameter(30),parameter(31),parameter(32),parameter(33));

     %pick start point of window
     start1=20;
     start2=10;
     start3=240;
     %save concentration for 20 min windows 
     phase1=[phase1;concentration1(:,start1:start1+19)];
     phase2=[phase2;concentration2(:,start2:start2+19)];
     phase3=[phase3;concentration2(:,start3:start3+19)];
end

%get dx
split_point=list_param*list_state;

phase1_up=phase1(1:split_point,:);
phase1_down=phase1(split_point+1:end,:);
dX_phase1=phase1_up-phase1_down;

phase2_up=phase2(1:split_point,:);
phase2_down=phase2(split_point+1:end,:);
dX_phase2=phase2_up-phase2_down;

phase3_up=phase3(1:split_point,:);
phase3_down=phase3(split_point+1:end,:);
dX_phase3=phase3_up-phase3_down;

%get dX/dP 
%get xi 
base_conc1=save_conc(:,start1:start1+19);
base_conc2=save_conc(:,start2:start2+19);
base_conc3=save_conc(:,start3:start3+19);

%set up collection array
DX_phase1=[];
DX_phase2=[];
DX_phase3=[];

state_index=1;

for index=1:list_param
    
    end_count=state_index+list_state-1;
    
    dXdP_phase1=(dX_phase1(state_index:end_count,:)/(save_param(index)*(delta1-delta2)));
    dXdP_phase2=(dX_phase2(state_index:end_count,:)/(save_param(index)*(delta1-delta2))).*(save_param(index)./(base_conc2(:,:)));
    dXdP_phase3=(dX_phase3(state_index:end_count,:)/(save_param(index)*(delta1-delta2))).*(save_param(index)./(base_conc3(:,:)));
    
    DX_phase1=[DX_phase1;dXdP_phase1];
    DX_phase2=[DX_phase2;dXdP_phase2];
    DX_phase3=[DX_phase3;dXdP_phase3];
    
    [state_index,list_datapoint]=size(dXdP_phase1);
    state_index=state_index+1;
end

%rearrange to get scaled sensitivity coefficients for all data points
sens_array1=[];
sens_array2=[];
sens_array3=[];

for time_index=1:list_datapoint
    
    state_index=1;
    sens_coeff1=[];
    sens_coeff2=[];
    sens_coeff3=[];
    
    for param_index=1:list_param
              
        for index=1:list_state
           
            end_count=state_index+list_state-1;

            temp1=DX_phase1(state_index:end_count,time_index);
            temp2=DX_phase2(state_index:end_count,time_index);
            temp3=DX_phase3(state_index:end_count,time_index);
        end
        
        sens_coeff1=[sens_coeff1 temp1];
        sens_coeff2=[sens_coeff2 temp2];     
        sens_coeff3=[sens_coeff3 temp3];

        state_index=end_count+1;
    end
    sens_array1=[sens_array1;sens_coeff1];
    sens_array2=[sens_array2;sens_coeff2];
    sens_array3=[sens_array3;sens_coeff3];
end

%the sensitivity coefficient arrays can now be split into coefficient
%arrays at particular data-points if required but this is not necessary for
%the next steps.

%--------------------------------------------------------------------------%
%--------------------------------------------------------------------------%
%%part(c)
%--------------------------------------------------------------------------%

%integrate s to obtain time-averages sensitivity array S
time_step=1.0;
trap1=0;
trap2=0;
trap3=0;

start_index1=1;

for time_index=1:list_datapoint-1 
    
    end_count=start_index1+list_state-1;
    discrete1_1=sens_array1(start_index1:end_count,:);
    discrete1_2=sens_array2(start_index1:end_count,:);
    discrete1_3=sens_array3(start_index1:end_count,:);
    
    start_index2=end_count+1;
    end_count=start_index2+list_state-1;
    discrete2_1=sens_array1(start_index2:end_count,:);
    discrete2_2=sens_array2(start_index2:end_count,:);
    discrete2_3=sens_array3(start_index2:end_count,:);
    
    trap1=trap1+0.5*(discrete1_1+discrete2_1)*time_step;
    trap2=trap2+0.5*(discrete1_2+discrete2_2)*time_step;
    trap3=trap3+0.5*(discrete1_3+discrete2_3)*time_step;
    
    start_index1=start_index2;
end 

window_length=20;

S1=trap1/window_length;
S2=trap2/window_length;
S3=trap3/window_length;

%Singular Value Decomposition
[U1,N1,V1]=svd(S1,'econ');
[U2,N2,V2]=svd(S2,'econ');
[U3,N3,V3]=svd(S3,'econ');

rank_phase(:,1)=abs(U1(:,1));
rank_phase(:,2)=abs(U2(:,1));
rank_phase(:,3)=abs(U3(:,1));

title=["species" "rank Phase1" "rank Phase2" "rank Phase2"];
species=["mRNA1";"mRNA2";"mRNA3";"protein1";"protein2";"protein3"];
result=[species rank_phase];
result=[title;result];
save('rank.mat','result','-mat');

%--------------------------------------------------------------------------%
%---------------------------end--------------------------------------------%