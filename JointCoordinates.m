function [BA_f,BA_r,BA_a,AP_f,AP_r,AP_a,BS_f,BS_r,BS_a,new_time, R_5th_M_z,frames] = JointCoordinates(dynamic_trial,static_trial,fix)

[sR5M,sRGT,sRLE,sRLO,sRLS,sT1,sRDS,sCentroid,time,sRME,sRMS,sRTR,sRCR,sR2M,sACB,sRAC] = JCSextract(static_trial,fix);

load([dynamic_trial '.mat']);

sR5M

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

[R_5th_M_x R_5th_M_y R_5th_M_z] = extract_XYZ(trial.pos_data.R5M);
[RLE_x RLE_y RLE_z] = extract_XYZ(trial.pos_data.RLE);
[R_Tricep_x R_Tricep_y R_Tricep_z] = extract_XYZ(trial.pos_data.RTR);
[ME_x ME_y ME_z] = extract_XYZ(trial.pos_data.RME);
[MS_x MS_y MS_z] = extract_XYZ(trial.pos_data.RMS);
[RGT_x, RGT_y, RGT_z] = extract_XYZ(trial.pos_data.RGT);
[RLO_x, RLO_y, RLO_z] = extract_XYZ(trial.pos_data.RLO);
[RLS_x, RLS_y, RLS_z] = extract_XYZ(trial.pos_data.RLS);
[T1_x T1_y T1_z] = extract_XYZ(trial.pos_data.T1);
[Centroid_x Centroid_y Centroid_z] = extract_XYZ(trial.pos_data.Centroid);
[RDS_x RDS_y RDS_z] = extract_XYZ(trial.pos_data.RDS);
[RCR_x RCR_y RCR_z] = extract_XYZ(trial.pos_data.RCR);
[RAC_x RAC_y RAC_z] = extract_XYZ(trial.pos_data.RAC);
[ACB_x ACB_y ACB_z] = extract_XYZ(trial.pos_data.ACB);


end_frame = length(trial.pos_data.R5M);

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


for k = 1:size(R_5th_M_x,2)
    
    RME_x(:,k) = R_Tricep_x(:,k) - static_Tri_x + static_ME_x;
    RME_y(:,k) = R_Tricep_y(:,k) - static_Tri_y + static_ME_y;
    RME_z(:,k) = R_Tricep_z(:,k) - static_Tri_z + static_ME_z;
    
    
    RMS_x(:,k) = RCR_x(:,k) - static_RCR_x + static_MS_x;
    RMS_y(:,k) = RCR_y(:,k) - static_RCR_y + static_MS_y;
    RMS_z(:,k) = RCR_z(:,k) - static_RCR_z + static_MS_z;
    
    R2M_x(:,k) = R_5th_M_x(:,k) - static_5th_M_x + static_R_2_x;
    R2M_y(:,k) = R_5th_M_y(:,k) - static_5th_M_y + static_R_2_y;
    R2M_z(:,k) = R_5th_M_z(:,k) - static_5th_M_z + static_R_2_z;
    
end

%%

%The following code uses an unweighted least squares function to
%approximate the translation vector r and rotation matrix R.

static_RME = [static_ME_x(1) static_ME_y(1) static_ME_z(1)]; %x-y-z coordinates of Medial Epicondyle
                                                   %at time t1=0 (static)

RGT_1 = [sRGT(1,1) sRGT(1,2) sRGT(1,3)]; %x-y-z coordinates of Greater Tubercle
                                                   %at time t1=0 (static)
    
RLE_1 = [sRLE(1,1) sRLE(1,2) sRLE(1,3)]; %x-y-z coordinates of Lateral Epicondyle
                                                   %at time t1=0 (static)
    
R_Tricep_1 = [sRTR(1,1) sRTR(1,2) sRTR(1,3)]; %x-y-z coordinates of Tricep
                                                   %at time t1=0 (static)
m = 3; %Number of markers 

center_of_rot_1 = 1/m*(RGT_1 + RLE_1 + R_Tricep_1); %Calculates the center
%of rotation at time t1=0 (static) for the following markers: Greater
%Tubercle, Lateral Epicondyle and Tricep.

A = 1/m*((RGT_1 - center_of_rot_1)'*(RGT_1 - center_of_rot_1) +... Calculates the distribution
        (RLE_1 - center_of_rot_1)'*(RLE_1 - center_of_rot_1) +... %matrix A at time t1=0
        (R_Tricep_1 - center_of_rot_1)'*(R_Tricep_1 - center_of_rot_1));

for j = 1:size(RME_x,2)
    
    for i = 1:length(RME_x)
    
    m = 3; %Number of markers    
        
    RGT_2 = [RGT_x(i) RGT_y(i) RGT_z(i)]; %x-y-z coordinates of Greater Tubercle
                                                   %at time t2 (dynamic)
    
    RLE_2 = [RLE_x(i) RLE_y(i) RLE_z(i)]; %x-y-z coordinates of Lateral Epicondyle
                                                   %at time t2 (dynamic)
    
    R_Tricep_2 = [R_Tricep_x(i) R_Tricep_y(i) R_Tricep_z(i)]; %x-y-z coordinates of tricep
                                                   %at time t2 (dynamic)
        
    center_of_rot_2 = 1/m*(RGT_2 + RLE_2 + R_Tricep_2); %Calculates the center
%of rotation at time t2 (dynamic) for the following markers: Greater
%Tubercle, Lateral Epicondyle and Tricep.
    
    G = 1/m*((RGT_2 - center_of_rot_2)'*(RGT_1 - center_of_rot_1)... %Calculates the distribution
        + (RLE_2 - center_of_rot_2)'*(RLE_1 - center_of_rot_1)... %matrix G at time t2 (dynamic)
        + (R_Tricep_2 - center_of_rot_2)'*(R_Tricep_1 - center_of_rot_1));
    
    r = center_of_rot_2 - center_of_rot_1; %Calculates translation vector from time t1 to t2
    
    R = G*inv(A); %Calculates rotation matrix from time t1 to t2
    
    RME_2 = [RME_x(i) RME_y(i) RME_z(i)]; %x-y-z coordinates of Medial Epicondyle
                                          %at time t2 (dynamic)
    
    RME = RME_2' + r' + R*(RME_2' - static_RME'); %Calculates the approximated position
    %of the Medial Epicondyle
    
    new_RME_x(i) = RME(1); %saves RME x data to a column vector for each iteration
    
    new_RME_y(i) = RME(2); %saves RME y data to a column vector for each iteration
    
    new_RME_z(i) = RME(3); %saves RME z data to a column vector for each iteration
    
    end
    
end

new_time = time(1:length(RGT_x));

% figure(1)
% plot(new_time,RME_z(:,1),'r')
% hold on
% plot(new_time,RLE_z(:,1),'b')
% hold off
% xlabel('time(sec)')
% ylabel('position(m)')
% title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (z-axis)')
% legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% % axis([0 0.5 0 0.35])
% 
% figure(2)
% plot(new_time,RME_x(:,1),'r')
% hold on
% plot(new_time,RLE_x(:,1),'b')
% hold off
% xlabel('time(sec)')
% ylabel('position(m)')
% title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (x-axis)')
% legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% % axis([0 0.5 0 0.35])
% 
% figure(3)
% plot(new_time,RME_y(:,1),'r')
% hold on
% plot(new_time,RLE_y(:,1),'b')
% hold off
% xlabel('time(sec)')
% ylabel('position(m)')
% title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (y-axis)')
% legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% % axis([0 0.5 0 0.35])


%%

g = 1;

v_static = [(static_ME_x - static_RLE_x) (static_ME_y - static_RLE_y) (static_ME_z - static_RLE_z)];

d_static = sqrt(v_static(:,1).^2 + v_static(:,2).^2 + v_static(:,3).^2);

v = [(RME_x(:,g) - RLE_x(:,g)) (RME_y(:,g) - RLE_y(:,g)) (RME_z(:,g) - RLE_z(:,g))];
vme = [(RME_x(:,g) ) (RME_y(:,g) ) (RME_z(:,g) )];
vle = [( RLE_x(:,g)) (RLE_y(:,g)) ( RLE_z(:,g))];

d = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
d = d/10;
dme = sqrt(vme(:,1).^2 + vme(:,2).^2 + vme(:,3).^2);
dle = sqrt(vle(:,1).^2 + vle(:,2).^2 + vle(:,3).^2);

% 
% figure(4)
% plot(d_static,'r')    
% hold on
% plot(d,'b')
% hold off
% xlabel('Time(sec)')
% ylabel('Distance(m)')
% title('Distance Between Medial and Lateral Epicondyles')
% legend({'Static Length','Dynamic Length'},'location','southeast')
% 
% figure(122)
% plot(new_time,dme,'r')
% hold on
% plot(new_time,dle,'b')
% hold off
% xlabel('Time(sec)')
% ylabel('mag(m)')
% title(' Medial and Lateral Epicondyles')
% % axis([0 0.5 0 0.2])
% legend({'me','le'},'location','southeast')
% 

%Establish Joint Coordinate System for BRACHIUM and ANTIBRACHIUM
for b = 1:size(RGT_x,2)
    
    %Local Coordinate System for PROXIMAL segment
  
    %z-axis
    
    z_prox_num = [(RLE_x(:,b) - RME_x(:,b)) (RLE_y(:,b) - RME_y(:,b)) (RLE_z(:,b) - RME_z(:,b))];
    
    z_prox_denom = sqrt(z_prox_num(:,1).^2 + z_prox_num(:,2).^2 + z_prox_num(:,3).^2);
    
    z_prox = z_prox_num./z_prox_denom;
    
    %x-axis
    
    RGT_RLE = [(RGT_x(:,b) - RLE_x(:,b)) (RGT_y(:,b) - RLE_y(:,b)) (RGT_z(:,b) - RLE_z(:,b))];
    
    x_prox_num = cross(RGT_RLE,z_prox);
    
    x_prox_denom = sqrt(x_prox_num(:,1).^2 + x_prox_num(:,2).^2 + x_prox_num(:,3).^2);
    
    x_prox = x_prox_num./x_prox_denom;
    
    %y-axis
    
    y_prox = cross(z_prox,x_prox);
    
    
    %Local Coordinate System for DISTAL segment
    
    %z-axis
    
    z_dist_num = [(RLS_x(:,b) - RMS_x(:,b)) (RLS_y(:,b) - RMS_y(:,b)) (RLS_z(:,b) - RMS_z(:,b))];
    
    z_dist_denom = sqrt(z_dist_num(:,1).^2 + z_dist_num(:,2).^2 + z_dist_num(:,3).^2);
    
    z_dist = z_dist_num./z_dist_denom;
    
    %x-axis
    
    RLO_RLS = [(RLO_x(:,b) - RLS_x(:,b)) (RLO_y(:,b) - RLS_y(:,b)) (RLO_z(:,b) - RLS_z(:,b))];
    
    x_dist_num = cross(RLO_RLS,z_dist);
    
    x_dist_denom = sqrt(x_dist_num(:,1).^2 + x_dist_num(:,2).^2 + x_dist_num(:,3).^2);
    
    x_dist = x_dist_num./x_dist_denom;
    
    %y-axis
    
    y_dist = cross(z_dist,x_dist);

    %Setting up the Joint Coordinate System
    
    FA = cross(y_dist,z_prox); %Calculates the floating axis
    
    alpha(:,b) = acosd(dot(x_prox,FA,2)); %Calculates the flexion angle in degrees
    
    gamma(:,b) = acosd(dot(x_dist,FA,2)); %Calculates internal/external rotation in degrees
    
    beta(:,b) = acosd(dot(z_prox,y_dist,2)); %Calculates abduction angle in degrees
    
    %///   
end

%Establish Joint Coordinate System for the ANTIBRACHIUM and PAW
%Local Coordinate System for PROXIMAL segment
    
    %z-axis
    
    z_prox_2 = z_dist; 
    
    %x-axis
    
    x_prox_2 = x_dist;
    
    %y-axis
    
    y_prox_2 = y_dist;
    
    %Local Coordinate System for Distal segment
    
    %z-axis
    
    z_dist_2_num = [(R_5th_M_x(:,b) - R2M_x(:,b)) (R_5th_M_y(:,b) - R2M_y(:,b)) (R_5th_M_z(:,b) - R2M_z(:,b))];
    
    z_dist_2_denom = sqrt(z_dist_2_num(:,1).^2 + z_dist_2_num(:,2).^2 + z_dist_2_num(:,3).^2);
    
    z_dist_2 = z_dist_2_num./z_dist_2_denom;
    
    %x-axis
    
    RAC_R_5th_M = [(ACB_x(:,b) - R_5th_M_x(:,b)) (ACB_y(:,b) - R_5th_M_y(:,b)) (ACB_z(:,b) - R_5th_M_z(:,b))];
    
    x_dist_2_num = cross(RAC_R_5th_M,z_dist_2);
    
    x_dist_2_denom = sqrt(x_dist_2_num(:,1).^2 + x_dist_2_num(:,2).^2 + x_dist_2_num(:,3).^2);
    
    x_dist_2 = x_dist_2_num./x_dist_2_denom;
    
    %y-axis
    
    y_dist_2 = cross(z_dist_2,x_dist_2);
    
    %Floating Axis
    
    FA_2 = cross(y_dist_2,z_prox_2);
    
    %Range of Motion of Joint
    
    alpha_2(:,b) = acosd(dot(x_prox_2,FA_2,2)); %Calculates the flexion angle in degrees
    
    gamma_2(:,b) = acosd(dot(x_dist_2,FA_2,2)); %Calculates internal/external rotation in degrees
    
    beta_2(:,b) = acosd(dot(z_prox_2,y_dist_2,2)); %Calculates abduction angle in degrees
    
    
% Not all segments have the same length for some reason which is causing
% prox_x from the BA to have an irregular length
%Establish Joint Coordinate System for BRACHIUM and SCAPULA
    
%Local Coordinate System for PROXIMAL segment
  
%z-axis

z_prox_3_num = [(RAC_x(:,b) - T1_x(:,b)) (RAC_y(:,b) - T1_y(:,b)) (RAC_z(:,b) - T1_z(:,b))];
    
z_prox_3_denom = sqrt(z_prox_3_num(:,1).^2 + z_prox_3_num(:,2).^2 + z_prox_3_num(:,3).^2);
    
z_prox_3 = z_prox_3_num./z_prox_3_denom;
    
%x-axis
    
R_Acrom_RDS = [(RDS_x(:,b) - T1_x(:,b)) (RDS_y(:,b) - T1_y(:,b)) (RDS_z(:,b) - T1_z(:,b))];
    
    x_prox_3_num = cross(R_Acrom_RDS,z_prox_3);
    
    x_prox_3_denom = sqrt(x_prox_3_num(:,1).^2 + x_prox_3_num(:,2).^2 + x_prox_3_num(:,3).^2);
    
    x_prox_3 = x_prox_3_num./x_prox_3_denom;
    
    %y-axis
    
    y_prox_3 = cross(z_prox_3,x_prox_3);

    %Local Coordinate System for DISTAL segment
    
    %z-axis
    
    z_dist_3 = z_prox;
    
    %x-axis
    
    x_dist_3 = x_prox;
    
    %y-axis
    
    y_dist_3 = y_prox;

    
    %Floating Axis
    
    FA_3 = cross(y_dist_3,z_prox_3);
    
    %Range of Motion of joint
    
    alpha_3(:,b) = acosd(dot(x_prox_3,FA_3,2)); %Calculates the flexion angle in degrees
    
    gamma_3(:,b) = acosd(dot(x_dist_3,FA,2)); %Calculates internal/external rotation in degrees
    
    beta_3(:,b) = acosd(dot(z_prox_3,y_dist_3,2)); %Calculates abduction angle in degrees
    

g2 = 1; 

BA_f = alpha(:,g2);
AP_f = alpha_2(:,g2);
BS_f = alpha_3(:,g2);

BA_r = gamma(:,g2);
AP_r = gamma_2(:,g2);
BS_r = gamma_3(:,g2);

BA_a = beta(:,g2);
AP_a = beta_2(:,g2);
BS_a = beta_3(:,g2);



%Graphs the velocity of T1 during each gait cycle

s = trial.gait_cycles(1).start_frame;
e = trial.gait_cycles(1).end_frame;

vT1 = T1_z(s:e);
vel = [];


for i=1:length(vT1)-1
    vel(i) = abs(((vT1(i+1) - vT1(i)) / 0.01) / 1000);
end

figure(9)
plot((s:e-1)/100,vel)
xlabel('seconds(s)')
ylabel('velocity(m/s)')
title('T1_z velocity during gc#1')

gait_cycle_count

for i = 1:gait_cycle_count
    
    if(~trial.gait_cycles(i).is_good)
        continue;
    end
    
    
    start_frame = trial.gait_cycles(i).start_frame;
    end_frame = trial.gait_cycles(i).end_frame;
    
    tBA_f = BA_f(start_frame:end_frame);
    tBA_r = BA_r(start_frame:end_frame);
    tBA_a = BA_a(start_frame:end_frame);
    
    tAP_f = AP_f(start_frame:end_frame);
    tAP_r = AP_r(start_frame:end_frame);
    tAP_a = AP_a(start_frame:end_frame);
    
    tBS_f = BS_f(start_frame:end_frame);
    tBS_r = BS_r(start_frame:end_frame);
    tBS_a = BS_a(start_frame:end_frame);
    
   
    frames(i) = end_frame - start_frame;
    x(i) = struct("ba_f",tBA_f,"ba_r",tBA_r,"ba_a",tBA_a,"ap_f",tAP_f,"ap_r",tAP_r,"ap_a",tAP_a,"bs_f",tBS_f,"bs_r",tBS_r,"bs_a",tBS_a,"start_frame",start_frame,"end_frame",end_frame);

end

jc_angles = struct("trial_name",trial.trial_name,"angles",x,"ba_f",BA_f,"ba_r",BA_r,"ba_a",BA_a,"ap_f",AP_f,"ap_r",AP_r,"ap_a",AP_a,"bs_f",BS_f,"bs_r",BS_r,"bs_a",BS_a);

temp_name = char(trial.trial_name + " Angles");

save([temp_name '.mat'],'jc_angles')

R_5th_M_z = R_5th_M_z * 1000;


%Handles exporting the angle data to spreadsheets

t = struct2table(jc_angles.angles,'AsArray',true);
y = table2array(t);
y_1 = y;

writecell(y_1,'chester_angles.xlsx');

% ba_table = [];
% 
% for i = 1:gait_cycle_count 
%     
%     if(~isempty(jc_angles.angles(i).ba_f))
%         ba_table = table(
%     end
% end



        
    








end








