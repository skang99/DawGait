
%Load Gait Cycles
load 'Wesson Trial 4 No Vest.mat'

%Read in static data

[static_data,static_COOR] = xlsread('Trials/Canin Static Trial No Vest.xlsx');

%%

%The following loop creates matrices of the x, y and z position data of the
%Medial Epicondyle during the dynamics trials using their respective static data.

for k = 1:size(RGT_x_GC,2)
    
    RME_x(:,k) = R_Tricep_x_GC(:,k) - static_Tri_x + static_ME_x;
    RME_y(:,k) = R_Tricep_y_GC(:,k) - static_Tri_y + static_ME_y;
    RME_z(:,k) = R_Tricep_z_GC(:,k) - static_Tri_z + static_ME_z;
    
    
    RMS_x(:,k) = RCR_x_GC(:,k) - static_RCR_x + static_MS_x;
    
    RMS_y(:,k) = RCR_y_GC(:,k) - static_RCR_y + static_MS_y;
    
    RMS_z(:,k) = RCR_z_GC(:,k) - static_RCR_z + static_MS_z;
    
    R_2nd_M_x(:,k) = R_5th_M_x_GC(:,k) - static_5th_M_x + static_2nd_M_x;
    
    R_2nd_M_y(:,k) = R_5th_M_y_GC(:,k) - static_5th_M_y + static_2nd_M_y;
    
    R_2nd_M_z(:,k) = R_5th_M_z_GC(:,k) - static_5th_M_z + static_2nd_M_z;
    
end


%The following code uses an unweighted least squares function to
%approximate the translation vector r and rotation matrix R.

static_RME = [static_ME_x(1) static_ME_y(1) static_ME_z(1)]; %x-y-z coordinates of Medial Epicondyle
                                                   %at time t1=0 (static)

RGT_1 = [static_RGT(1,1) static_RGT(1,2) static_RGT(1,3)]; %x-y-z coordinates of Greater Tubercle
                                                   %at time t1=0 (static)
    
RLE_1 = [static_RLE(1,1) static_RLE(1,2) static_RLE(1,3)]; %x-y-z coordinates of Lateral Epicondyle
                                                   %at time t1=0 (static)
    
R_Tricep_1 = [static_Tricep(1,1) static_Tricep(1,2) static_Tricep(1,3)]; %x-y-z coordinates of Tricep
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
        
    RGT_2 = [RGT_x_GC(i) RGT_y_GC(i) RGT_z_GC(i)]; %x-y-z coordinates of Greater Tubercle
                                                   %at time t2 (dynamic)
    
    RLE_2 = [RLE_x_GC(i) RLE_y_GC(i) RLE_z_GC(i)]; %x-y-z coordinates of Lateral Epicondyle
                                                   %at time t2 (dynamic)
    
    R_Tricep_2 = [R_Tricep_x_GC(i) R_Tricep_y_GC(i) R_Tricep_z_GC(i)]; %x-y-z coordinates of tricep
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


new_time = time(1:length(RGT_x_GC))

figure(1)
plot(new_time,RME_z(:,1),'r')
hold on
plot(new_time,RLE_z_GC(:,1),'b')
hold off
xlabel('time(sec)')
ylabel('position(m)')
title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (z-axis)')
legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% axis([0 0.5 0 0.35])

figure(2)
plot(new_time,RME_x(:,1),'r')
hold on
plot(new_time,RLE_x_GC(:,1),'b')
hold off
xlabel('time(sec)')
ylabel('position(m)')
title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (x-axis)')
legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% axis([0 0.5 0 0.35])

figure(3)
plot(new_time,RME_y(:,1),'r')
hold on
plot(new_time,RLE_y_GC(:,1),'b')
hold off
xlabel('time(sec)')
ylabel('position(m)')
title('Dynamic Motion of Medial Epicondyle vs. Lateral Epicondyle (y-axis)')
legend({'Medial Epicondyle','Lateral Epicondyle'},'location','southeast')
% axis([0 0.5 0 0.35])

%%

g = 1;

v_static = [(static_ME_x - static_RLE_x) (static_ME_y - static_RLE_y) (static_ME_z - static_RLE_z)];

d_static = sqrt(v_static(:,1).^2 + v_static(:,2).^2 + v_static(:,3).^2);

v = [(RME_x(:,g) - RLE_x_GC(:,g)) (RME_y(:,g) - RLE_y_GC(:,g)) (RME_z(:,g) - RLE_z_GC(:,g))];
vme = [(RME_x(:,g) ) (RME_y(:,g) ) (RME_z(:,g) )];
vle = [( RLE_x_GC(:,g)) (RLE_y_GC(:,g)) ( RLE_z_GC(:,g))];

d = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
dme = sqrt(vme(:,1).^2 + vme(:,2).^2 + vme(:,3).^2);
dle = sqrt(vle(:,1).^2 + vle(:,2).^2 + vle(:,3).^2);
figure(4)
plot(new_time,d_static,'r')
hold on
plot(new_time,d,'b')
hold off
xlabel('Time(sec)')
ylabel('Distance(m)')
title('Distance Between Medial and Lateral Epicondyles')
 axis([0 0.5 0 0.2])
legend({'Static Length','Dynamic Length'},'location','southeast')

figure(122)
plot(new_time,dme,'r')
hold on
plot(new_time,dle,'b')
hold off
xlabel('Time(sec)')
ylabel('mag(m)')
title(' Medial and Lateral Epicondyles')
% axis([0 0.5 0 0.2])
legend({'me','le'},'location','southeast')

%%

%Establish Joint Coordinate System for BRACHIUM and ANTIBRACHIUM


for b = 1:size(RGT_x_GC,2)
    
    %Local Coordinate System for PROXIMAL segment
  
    %z-axis
    
    z_prox_num = [(RLE_x_GC(:,b) - RME_x(:,b)) (RLE_y_GC(:,b) - RME_y(:,b)) (RLE_z_GC(:,b) - RME_z(:,b))];
    
    z_prox_denom = sqrt(z_prox_num(:,1).^2 + z_prox_num(:,2).^2 + z_prox_num(:,3).^2);
    
    z_prox = z_prox_num./z_prox_denom;
    
    %x-axis
    
    RGT_RLE = [(RGT_x_GC(:,b) - RLE_x_GC(:,b)) (RGT_y_GC(:,b) - RLE_y_GC(:,b)) (RGT_z_GC(:,b) - RLE_z_GC(:,b))];
    
    x_prox_num = cross(RGT_RLE,z_prox);
    
    x_prox_denom = sqrt(x_prox_num(:,1).^2 + x_prox_num(:,2).^2 + x_prox_num(:,3).^2);
    
    x_prox = x_prox_num./x_prox_denom;
    
    %y-axis
    
    y_prox = cross(z_prox,x_prox);
    
    
    %Local Coordinate System for DISTAL segment
    
    %z-axis
    
    z_dist_num = [(RLS_x_GC(:,b) - RMS_x(:,b)) (RLS_y_GC(:,b) - RMS_y(:,b)) (RLS_z_GC(:,b) - RMS_z(:,b))];
    
    z_dist_denom = sqrt(z_dist_num(:,1).^2 + z_dist_num(:,2).^2 + z_dist_num(:,3).^2);
    
    z_dist = z_dist_num./z_dist_denom;
    
    %x-axis
    
    RLO_RLS = [(RLO_x_GC(:,b) - RLS_x_GC(:,b)) (RLO_y_GC(:,b) - RLS_y_GC(:,b)) (RLO_z_GC(:,b) - RLS_z_GC(:,b))];
    
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
    
    z_dist_2_num = [(R_5th_M_x_GC(:,b) - R_2nd_M_x(:,b)) (R_5th_M_y_GC(:,b) - R_2nd_M_y(:,b)) (R_5th_M_z_GC(:,b) - R_2nd_M_z(:,b))];
    
    z_dist_2_denom = sqrt(z_dist_2_num(:,1).^2 + z_dist_2_num(:,2).^2 + z_dist_2_num(:,3).^2);
    
    z_dist_2 = z_dist_2_num./z_dist_2_denom;
    
    %x-axis
    
    RAC_R_5th_M = [(RAC_x_GC(:,b) - R_5th_M_x_GC(:,b)) (RAC_y_GC(:,b) - R_5th_M_y_GC(:,b)) (RAC_z_GC(:,b) - R_5th_M_z_GC(:,b))];
    
    x_dist_2_num = cross(RAC_R_5th_M,z_dist_2);
    
    x_dist_2_denom = sqrt(x_dist_2_num(:,1).^2 + x_dist_2_num(:,2).^2 + x_dist_2_num(:,3).^2);
    
    x_dist_2 = x_dist_2_num./x_dist_2_denom;
    
    %y-axis
    
    y_dist_2 = cross(z_dist_2,x_dist_2);
    
    
    %Floating Access
    
    FA_2 = cross(y_dist_2,z_prox_2);
    
    %Range of Motion of Joint
    
    alpha_2(:,b) = acosd(dot(x_prox_2,FA_2,2)); %Calculates the flexion angle in degrees
    
    gamma_2(:,b) = acosd(dot(x_dist_2,FA_2,2)); %Calculates internal/external rotation in degrees
    
    beta_2(:,b) = acosd(dot(z_prox_2,y_dist_2,2)); %Calculates abduction angle in degrees
    
    
    
    %Establish Joint Coordinate System for BRACHIUM and SCAPULA
    
    %Local Coordinate System for PROXIMAL segment
  
    %z-axis
    
    z_prox_3_num = [(R_Acrom_x_GC(:,b) - L_Acrom_x_GC(:,b)) (R_Acrom_y_GC(:,b) - L_Acrom_y_GC(:,b)) (R_Acrom_z_GC(:,b) - L_Acrom_z_GC(:,b))];
    
    z_prox_3_denom = sqrt(z_prox_3_num(:,1).^2 + z_prox_3_num(:,2).^2 + z_prox_3_num(:,3).^2);
    
    z_prox_3 = z_prox_3_num./z_prox_3_denom;
    
    
    %x-axis
    
    R_Acrom_RDS = [(RDS_x_GC(:,b) - R_Acrom_x_GC(:,b)) (RDS_y_GC(:,b) - R_Acrom_y_GC(:,b)) (RDS_z_GC(:,b) - R_Acrom_z_GC(:,b))];
    
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
   

    

g2 = 2;

figure(5)
plot(new_time,alpha(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Brachium/Antibrachium Flexion')
%axis([0 0.5 0 125])
legend({'Flexion Angle','Gait Cycle'},'location','southeast')


figure(6)
plot(new_time,alpha_2(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Antibrachium/Paw Flexion')
%axis([0 0.5 0 125])
legend({'Flexion Angle','Gait Cycle'},'location','southeast')


figure(7)
plot(new_time,alpha_3(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Scapula/Brachium Flexion')
%axis([0 0.5 0 125])
legend({'Flexion Angle','Gait Cycle'},'location','southeast')

figure(8)
plot(new_time,gamma(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Brachium/Antibrachium Int/Ext Rotation')
%axis([0 0.5 0 125])
legend({'Int Rotation Angle','Gait Cycle'},'location','southeast')


figure(9)
plot(new_time,gamma_2(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Antibrachium/Paw Int/Ext Rotation')
%axis([0 0.5 0 125])
legend({'Int Rotation Angle','Gait Cycle'},'location','southeast')


figure(10)
plot(new_time,gamma_3(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Scapula/Brachium Int/Ext Rotation')
%axis([0 0.5 0 125])
legend({'Int Rotation Angle','Gait Cycle'},'location','southeast')

figure(11)
plot(new_time,beta(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Brachium/Antibrachium Abduction/Adduction')
%axis([0 0.5 0 125])
legend({'Abduction Angle','Gait Cycle'},'location','southeast')


figure(12)
plot(new_time,beta_2(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Antibrachium/Paw Abduction/Adduction')
%axis([0 0.5 0 125])
legend({'Abduction Angle','Gait Cycle'},'location','southeast')


figure(13)
plot(new_time,beta_3(:,g2),'r')
hold on
plot(new_time,R_5th_M_z_GC(:,g2)*1000,'b')
hold off
xlabel('time(sec)')
ylabel('angle(deg)')
title('Scapula/Brachium Abduction/Adduction')
%axis([0 0.5 0 125])
legend({'Abduction Angle','Gait Cycle'},'location','southeast')



% figure(4)
%  for i = 1:length(new_RME_x)
% 
%      
%     plot([RGT_x_GC(i),RLE_x_GC(i)],[RGT_z_GC(i),RLE_z_GC(i)],'r-o','Linewidth',1.5)
%     
%     hold on
%     
%     plot([RLE_x_GC(i),R_Tricep_x_GC(i)],[RLE_z_GC(i),R_Tricep_z_GC(i)],'b-o','Linewidth',1.5)
%     
%     hold on
% 
%     plot([RGT_x_GC(i),R_Tricep_x_GC(i)],[RGT_z_GC(i),R_Tricep_z_GC(i)],'y-o','Linewidth',1.5)
%     
%     hold on
%     
%     plot([R_Tricep_x_GC(i),new_RME_x(i)],[R_Tricep_z_GC(i),new_RME_z(i)],'r-o','Linewidth',1.5)
%     
%     hold off
%     
%     axis([-1 2 0 0.5])
%     
%     text(RGT_x_GC(i),RGT_z_GC(i),'RGT')
%     text(RLE_x_GC(i),RLE_z_GC(i),'RLE')
%     text(R_Tricep_x_GC(i),R_Tricep_z_GC(i),'R Tri')
%     pause(0.01)
%     
%     
%  end





