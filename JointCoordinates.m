function [BA_f,BA_r,BA_a,AP_f,AP_r,AP_a,BS_f,BS_r,BS_a,new_time, R_5th_M_z,frames,no_good_cycles,gait_cycle_count] = JointCoordinates(dynamic_trial,static_trial)

try
[sR5M,sRGT,sRLE,sRLO,sRLS,sT1,sRDS,sCentroid,time,sRME,sRMS,sRTR,sRCR,sR2M,sACB,sRAC] = JCSextract(static_trial);
catch exception
    rethrow(exception)
end

load(['Produced Data/' dynamic_trial '.mat'])

[static_5th_M_x static_5th_M_y static_5th_M_z] = extract_XYZ(sR5M);
[static_RLE_x static_RLE_y static_RLE_z] = extract_XYZ(sRLE);
[static_Tri_x static_Tri_y static_Tri_z] = extract_XYZ(sRTR);
[static_ME_x static_ME_y static_ME_z] = extract_XYZ(sRME);
[static_RCR_x static_RCR_y static_RCR_z] = extract_XYZ(sRCR);
[static_MS_x static_MS_y static_MS_z] = extract_XYZ(sRMS);
[static_R_2_x static_R_2_y static_R_2_z] = extract_XYZ(sR2M);
[static_T1_x static_T1_y static_T1_z] = extract_XYZ(sT1);
[static_RAC_x static_RAC_y static_RAC_z] = extract_XYZ(sRAC);
[static_Centroid_x static_Centroid_y static_Centroid_z] = extract_XYZ(sCentroid);
[static_RDS_x static_RDS_y static_RDS_z] = extract_XYZ(sRDS);

gait_cycle_count = length(trial.gait_cycles);

no_good_cycles = 0;


R5M = trial.pos_data.R5M;
[R5M_x R5M_y R5M_z] = extract_XYZ(R5M);
RLE = trial.pos_data.RLE;
RTR = trial.pos_data.RTR;
RME = trial.pos_data.RME;
RMS = trial.pos_data.RMS;
RGT = trial.pos_data.RGT;
RLO = trial.pos_data.RLO;
RLS = trial.pos_data.RLS;
T1 = trial.pos_data.T1;
Centroid = trial.pos_data.Centroid;
RDS = trial.pos_data.RDS;
RCR = trial.pos_data.RCR;
RAC = trial.pos_data.RAC;
ACB = trial.pos_data.ACB;

%allows for static and dynamic matrix dimensions to agree
end_frame = length(trial.pos_data);

static_RLE_x = static_RLE_x(1:end_frame);
static_RLE_y = static_RLE_y(1:end_frame);
static_RLE_z = static_RLE_z(1:end_frame);

static_Tri_x = static_Tri_x(1:end_frame);
static_Tri_y = static_Tri_y(1:end_frame);
static_Tri_z = static_Tri_z(1:end_frame);

static_ME_x = static_ME_x(1:end_frame);
static_ME_y = static_ME_y(1:end_frame);
static_ME_z = static_ME_z(1:end_frame);

static_RCR_x = static_RCR_x(1:end_frame);
static_RCR_y = static_RCR_y(1:end_frame);
static_RCR_z = static_RCR_z(1:end_frame);

static_MS_x = static_MS_x(1:end_frame);
static_MS_y = static_MS_y(1:end_frame);
static_MS_z = static_MS_z(1:end_frame);

static_5th_M_x = static_5th_M_x(1:end_frame);
static_5th_M_y = static_5th_M_y(1:end_frame);
static_5th_M_z = static_5th_M_z(1:end_frame);


static_R_2_x = static_R_2_x(1:end_frame);
static_R_2_y = static_R_2_y(1:end_frame);
static_R_2_z = static_R_2_z(1:end_frame);


%%

%The following loop creates matrices of the x, y and z position data of the
%Medial Epicondyle during the dynamics trials using their respective static data.
RMS_x = [];
RMS_y = [];
RMS_z = [];

R2M_x = [];
R2M_y = [];
R2M_z = [];

[R_Tricep_x R_Tricep_y R_Tricep_z] = extract_XYZ(RTR);
[RCR_x RCR_y RCR_z] = extract_XYZ(RCR);

for k = 1:size(R5M_x,2)
    
    RME_x(:,k) = R_Tricep_x(:,k) - static_Tri_x + static_ME_x;
    RME_y(:,k) = R_Tricep_y(:,k) - static_Tri_y + static_ME_y;
    RME_z(:,k) = R_Tricep_z(:,k) - static_Tri_z + static_ME_z;
       
    RMS_x(:,k) = RCR_x(:,k) - static_RCR_x + static_MS_x;
    RMS_y(:,k) = RCR_y(:,k) - static_RCR_y + static_MS_y;
    RMS_z(:,k) = RCR_z(:,k) - static_RCR_z + static_MS_z;
    
    R2M_x(:,k) = R5M_x(:,k) - static_5th_M_x + static_R_2_x;
    R2M_y(:,k) = R5M_y(:,k) - static_5th_M_y + static_R_2_y;
    R2M_z(:,k) = R5M_z(:,k) - static_5th_M_z + static_R_2_z;
    
    
    
end

RME = [RME_x RME_y RME_z];
RMS = [RMS_x RMS_y RMS_z];
R2M = [R2M_x R2M_y RMS_z];

new_time = end_frame;

bm = [-1 -1 -1];

%Establish Joint Coordinate System for BRACHIUM and ANTIBRACHIUM
ps1 = 0;
ps2 = 0;
ps3 = 0;

[BA_f BA_r BA_a ps1 ps2 ps3] = create_angle_data(RLE,RME,RGT,RLS,RMS,RLO,bm,ps1,ps2,ps3);
[AP_f AP_r AP_a ps1 ps2 ps3] = create_angle_data(bm,bm,bm,R5M,R2M,ACB,bm,ps1,ps2,ps3);
[BS_f BS_r BS_a ps1 ps2 ps3] = create_angle_data(Centroid,T1,RDS,bm,bm,bm,bm,ps1,ps2,ps3);

% % BA_f = alpha(:,g2);
% % AP_f = alpha_2(:,g2);
% BS_f = alpha_3(:,g2);
% 
% % BA_r = gamma(:,g2);
% % AP_r = gamma_2(:,g2);
% BS_r = gamma_3(:,g2);
% 
% % BA_a = beta(:,g2);
% % AP_a = beta_2(:,g2);
% BS_a = beta_3(:,g2);

%Graphs the velocity of T1 during gc1

% s = trial.gait_cycles(1).start_frame;
% e = trial.gait_cycles(1).end_frame;

% vT1 = T1_y(s:e);
% vel = [];
% 
% 
% for i=1:length(vT1)-1
%     vel(i) = abs(((vT1(i+1) - vT1(i)) / 0.02) / 1000);
% end

% figure(9)
% plot((s:e-1)/100,vel)
% xlabel('seconds(s)')
% ylabel('velocity(m/s)')
% title('T1_y velocity during gc#1')

% %Show anyway? 
% if(no_good_cycles)
%     BA_f =0;BA_r=0;BA_a=0;AP_f=0;AP_r=0;AP_a=0;BS_f=0;BS_r=0;BS_a=0;new_time=0; R_5th_M_z=0;frames=0;
%     no_good_cycles = 1;
%     return;
% end

nBA_f = [];
nBA_r = [];
nBA_a = [];

nAP_f = [];
nAP_r = [];
nAP_a = [];

nBS_f = [];
nBS_r = [];
nBS_a = [];

for i = 1:gait_cycle_count
    
%     if(~trial.gait_cycles(i).is_good)
%         continue;
%     end
    
    
    start_frame = trial.gait_cycles(i).start_frame;
    end_frame = trial.gait_cycles(i).end_frame;
    gc_length = abs(start_frame - end_frame) + 1;
    
    if(trial.gait_cycles(i).RGT_RLE && trial.gait_cycles(i).RLO_RLS)
        %Elbow RGT_RLE & RLO_RLS
        disp("Creating elbow data for gc # " + i)
        tBA_f = BA_f(start_frame:end_frame);
        tBA_r = BA_r(start_frame:end_frame);
        tBA_a = BA_a(start_frame:end_frame);
        nBA_f = [nBA_f; tBA_f];
        nBA_r = [nBA_r; tBA_r];
        nBA_a = [nBA_a; tBA_a];
    else
        tBA_f = zeros(gc_length,1) + -1;
        tBA_r = zeros(gc_length,1) + -1;
        tBA_a = zeros(gc_length,1) + -1;
    end
    
    %Carpus RLO_RLS & ACB_MP5
    if(trial.gait_cycles(i).ACB_R5M && trial.gait_cycles(i).RLO_RLS)
        disp("Creating carpus data for gc # " + i)
        tAP_f = AP_f(start_frame:end_frame);
        tAP_r = AP_r(start_frame:end_frame);
        tAP_a = AP_a(start_frame:end_frame);
        nAP_f = [nAP_f; tAP_f];
        nAP_r = [nAP_r; tAP_r];
        nAP_a = [nAP_a; tAP_a];
    else
        tAP_f = zeros(gc_length,1) + -1;
        tAP_r = zeros(gc_length,1) + -1;
        tAP_a = zeros(gc_length,1) + -1;
    end 
    
    %shoulder
    if(trial.gait_cycles(i).RGT_RLE && trial.gait_cycles(i).RDS_Centroid && trial.gait_cycles(i).T1_Centroid)
        disp("Creating shoulder data for gc # " + i)
        tBS_f = BS_f(start_frame:end_frame); 
        tBS_r = BS_r(start_frame:end_frame);
        tBS_a = BS_a(start_frame:end_frame);
        nBS_f = [nBS_f; tBS_f];
        nBS_r = [nBS_r; tBS_r];
        nBS_a = [nBS_a; tBS_a];
    else
        tBS_f = zeros(gc_length,1) + -1;
        tBS_r = zeros(gc_length,1) + -1;
        tBS_a = zeros(gc_length,1) + -1;
    end
       
    %postcondition: all data included, excluding bad(seg check failure)/missing data
    x(i) = struct("ba_f",tBA_f,"ba_r",tBA_r,"ba_a",tBA_a,"ap_f",tAP_f,"ap_r",tAP_r,"ap_a",tAP_a,"bs_f",tBS_f,"bs_r",tBS_r,"bs_a",tBS_a,"start_frame",start_frame,"end_frame",end_frame);

end
frames = [];

for i = 1:length(trial.gait_cycles)
    frames(i) = trial.gait_cycles(i).start_frame;
    frames(i+1) = trial.gait_cycles(i).end_frame;
end



% nBS_f = nBS_f(~isnan(nBS_f));
% nBS_r = nBS_r(~isnan(nBS_r));
% nBS_a = nBS_a(~isnan(nBS_a));
% 
% nBA_f = nBA_f(~isnan(nBA_f));
% nBA_r = nBA_r(~isnan(nBA_r));
% nBA_a = nBA_a(~isnan(nBA_a));
% 
% nAP_f = nAP_f(~isnan(nAP_f));
% nAP_r = nAP_r(~isnan(nAP_r));
% nAP_a = nAP_a(~isnan(nAP_a));

jc_angles = struct("trial_name",trial.trial_name,"angles",x,"ba_f",nBA_f,"ba_r",nBA_r,"ba_a",nBA_a,"ap_f",nAP_f,"ap_r",nAP_r,"ap_a",nAP_a,"bs_f",nBS_f,"bs_r",nBS_r,"bs_a",nBS_a);

temp_name = char(trial.trial_name + " Angles");

save(['Produced Data/' temp_name '.mat'],'jc_angles')

R_5th_M_z = R5M_z;

name = {trial.trial_name
        ''
        'GC#'};

%TODO finish: create append
header = {'1' '' '' '' '' '' '' '' '' '2' '' '' '' '' '' '' '' '' '3' '' '' '' '' '' '' '' '' 
           'BA_f' 'BA_r' 'BA_a' 'AP_f' 'AP_r' 'AP_a' 'BS_f' 'BS_r' 'BS_a' 'BA_f' 'BA_r' 'BA_a' 'AP_f' 'AP_r' 'AP_a' 'BS_f' 'BS_r' 'BS_a' 'BA_f' 'BA_r' 'BA_a' 'AP_f' 'AP_r' 'AP_a' 'BS_f' 'BS_r' 'BS_a'};

filename = "Produced Data/" + trial.trial_name + " Angular Data.xlsx";

writecell(name,filename,'Range','A1')
writecell(header,filename,'Range','B3')

G1 = table(jc_angles.angles(1).ba_f,jc_angles.angles(1).bs_f,jc_angles.angles(1).ap_f,jc_angles.angles(1).ba_r,jc_angles.angles(1).bs_r,jc_angles.angles(1).ap_r,jc_angles.angles(1).ba_a,jc_angles.angles(1).bs_a,jc_angles.angles(1).ap_a);

if(gait_cycle_count > 1)
    G2 = table(jc_angles.angles(2).ba_f,jc_angles.angles(2).bs_f,jc_angles.angles(2).ap_f,jc_angles.angles(2).ba_r,jc_angles.angles(2).bs_r,jc_angles.angles(2).ap_r,jc_angles.angles(2).ba_a,jc_angles.angles(2).bs_a,jc_angles.angles(2).ap_a);
end

if(gait_cycle_count > 2)
    G3 = table(jc_angles.angles(3).ba_f,jc_angles.angles(3).bs_f,jc_angles.angles(3).ap_f,jc_angles.angles(3).ba_r,jc_angles.angles(3).bs_r,jc_angles.angles(3).ap_r,jc_angles.angles(3).ba_a,jc_angles.angles(3).bs_a,jc_angles.angles(3).ap_a);
end

writetable(G1,filename,'Sheet',1,'Range','B5','WriteVariableNames',false)

if(gait_cycle_count > 1)
    writetable(G2,filename,'Sheet',1,'Range','K5','WriteVariableNames',false)
end

if(gait_cycle_count > 2)
    writetable(G3,filename,'Sheet',1,'Range','T5','WriteVariableNames',false)
end

end








