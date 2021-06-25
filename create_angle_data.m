function [BA_f BA_r BA_a BS_f BS_r BS_a AP_f AP_r AP_a] = create_angle_data(RLE,RME,RGT,RLS,RMS,RLO,R5M,R2M,ACB,T1,RDS,RAC,Centroid)

[RLE_x RLE_y RLE_z] = extract_XYZ(RLE);
[RME_x RME_y RME_z] = extract_XYZ(RME);
[RGT_x RGT_y RGT_z] = extract_XYZ(RGT);
[RLS_x RLS_y RLS_z] = extract_XYZ(RLS);
[RMS_x RMS_y RMS_z] = extract_XYZ(RMS);
[RLO_x RLO_y RLO_z] = extract_XYZ(RLO);
[R5M_x R5M_y R5M_z] = extract_XYZ(R5M);
[R2M_x R2M_y R2M_z] = extract_XYZ(R2M);
[ACB_x ACB_y ACB_z] = extract_XYZ(ACB);
[T1_x T1_y T1_z] = extract_XYZ(T1);
[RDS_x RDS_y RDS_z] = extract_XYZ(RDS);
[RAC_x RAC_y RAC_z] = extract_XYZ(RAC);
[Centroid_x Centroid_y Centroid_z] = extract_XYZ(Centroid);


%Establish Joint Coordinate System for BRACHIUM and ANTIBRACHIUM or ELBOW

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
    
    FA = cross(z_prox,y_dist); %Calculates the floating axis
    
    alpha(:,b) = acosd(dot(FA,x_prox,2)); %Holds the flexion 
    gamma(:,b) = acosd(dot(FA,x_dist,2)); %Holds internal/external rotation
    beta(:,b) = acosd(dot(z_prox,y_dist,2)); %Holds abduction angle
     
end

%Establish Joint Coordinate System for the ANTIBRACHIUM and PAW or CARPUS
%Local Coordinate System for PROXIMAL segment
    
    %z-axis
    
    z_prox_2 = z_dist; 
    
    %x-axis
    
    x_prox_2 = x_dist;
    
    %y-axis
    
    y_prox_2 = y_dist;
    
    %Local Coordinate System for Distal segment
    
    %z-axis
    
    z_dist_2_num = [(R5M_x(:,b) - R2M_x(:,b)) (R5M_y(:,b) - R2M_y(:,b)) (R5M_z(:,b) - R2M_z(:,b))];
    
    z_dist_2_denom = sqrt(z_dist_2_num(:,1).^2 + z_dist_2_num(:,2).^2 + z_dist_2_num(:,3).^2);
    
    z_dist_2 = z_dist_2_num./z_dist_2_denom;
    
    %x-axis
    
    RAC_R5M = [(ACB_x(:,b) - R5M_x(:,b)) (ACB_y(:,b) - R5M_y(:,b)) (ACB_z(:,b) - R5M_z(:,b))];
    
    x_dist_2_num = cross(RAC_R5M,z_dist_2);
    
    x_dist_2_denom = sqrt(x_dist_2_num(:,1).^2 + x_dist_2_num(:,2).^2 + x_dist_2_num(:,3).^2);
    
    x_dist_2 = x_dist_2_num./x_dist_2_denom;
    
    %y-axis
    
    y_dist_2 = cross(z_dist_2,x_dist_2);
    
    %Floating Axis
    
    FA_2 = cross(y_dist_2,z_prox_2);
    
    %Range of Motion of Joint
    
    alpha_2(:,b) = asind(dot(-FA_2,y_prox_2,2)); %Holds the flexion 
    gamma_2(:,b) = acosd(dot(FA_2,x_dist_2,2)); %Holds internal/external rotation
    beta_2(:,b) = acosd(dot(z_prox_2,y_dist_2,2)); %Holds abduction angle
    

% Not all segments have the same length for some reason which is causing
% prox_x from the BA to have an irregular length

%Establish Joint Coordinate System for BRACHIUM and SCAPULA or SHOULDER
    
%Local Coordinate System for PROXIMAL segment
%What happened to using Centroid?
  
%z-axis

z_prox_3_num = [(T1_x(:,b)- Centroid_x(:,b)) (T1_y(:,b) - Centroid_y(:,b)) (T1_z(:,b) - Centroid_z(:,b))];
    
z_prox_3_denom = sqrt(z_prox_3_num(:,1).^2 + z_prox_3_num(:,2).^2 + z_prox_3_num(:,3).^2);
    
z_prox_3 = z_prox_3_num./z_prox_3_denom;
    
%x-axis
%Put back RDS_T1
    
Centroid_RDS = [(RDS_x(:,b) - Centroid_x(:,b)) (RDS_y(:,b) - Centroid_y(:,b)) (RDS_z(:,b) - Centroid_z(:,b))];
    
    x_prox_3_num = cross(Centroid_RDS,z_prox_3);
    
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
      
    alpha_3(:,b) = acosd(dot(FA_3,x_prox_3,2)); %Holds the flexion 
    gamma_3(:,b) = acosd(dot(FA_3,x_dist_3,2)); %Holds internal/external rotation
    beta_3(:,b) = acosd(dot(z_prox_3,y_dist_3,2)); %Holds abduction angle
    
    

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

end
