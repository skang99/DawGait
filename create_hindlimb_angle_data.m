function [BA_f BA_r BA_a BS_f BS_r BS_a AP_f AP_r AP_a] = create_hindlimb_angle_data(RLEP,RMEP,RGT,RLMA,RMMA,RFH,RMP5,R2M,RCAL,RISC,CRS,RIWG,PTC,side)

[RLEP_x RLEP_y RLEP_z] = extract_XYZ(RLEP);
[RMEP_x RMEP_y RMEP_z] = extract_XYZ(RMEP);
[RGT_x RGT_y RGT_z] = extract_XYZ(RGT);
[RLMA_x RLMA_y RLMA_z] = extract_XYZ(RLMA);
[RMMA_x RMMA_y RMMA_z] = extract_XYZ(RMMA);
[RFH_x RFH_y RFH_z] = extract_XYZ(RFH);
[RMP5_x RMP5_y RMP5_z] = extract_XYZ(RMP5);
[R2M_x R2M_y R2M_z] = extract_XYZ(R2M);
[RCAL_x RCAL_y RCAL_z] = extract_XYZ(RCAL);
[RISC_x RISC_y RISC_z] = extract_XYZ(RISC);
[CRS_x CRS_y CRS_z] = extract_XYZ(CRS);
[RIWG_x RIWG_y RIWG_z] = extract_XYZ(RIWG);
[PTC_x, PTC_y, PTC_z] = extract_XYZ(PTC);


for b = 1:size(RGT_x,2)   
    %Local Coordinate System for PROXIMAL segment
  
    %z-axis
    
    z_prox_num = [(RLEP_x(:,b) - RMEP_x(:,b)) (RLEP_y(:,b) - RMEP_y(:,b)) (RLEP_z(:,b) - RMEP_z(:,b))];
    
    z_prox_denom = sqrt(z_prox_num(:,1).^2 + z_prox_num(:,2).^2 + z_prox_num(:,3).^2);
    
    z_prox = z_prox_num./z_prox_denom;
    
    z_prox = z_prox * side;
    
    %x-axis
    
    RGT_RLE = [(RGT_x(:,b) - RLEP_x(:,b)) (RGT_y(:,b) - RLEP_y(:,b)) (RGT_z(:,b) - RLEP_z(:,b))];
    
    x_prox_num = cross(RGT_RLE,z_prox);
    
    x_prox_denom = sqrt(x_prox_num(:,1).^2 + x_prox_num(:,2).^2 + x_prox_num(:,3).^2);
    
    x_prox = x_prox_num./x_prox_denom;
    
    %y-axis
    
    y_prox = cross(z_prox,x_prox);
    
    
    %Local Coordinate System for DISTAL segment
    
    %z-axis
    
    z_dist_num = [(RLMA_x(:,b) - RMMA_x(:,b)) (RLMA_y(:,b) - RMMA_y(:,b)) (RLMA_z(:,b) - RMMA_z(:,b))];
    
    z_dist_denom = sqrt(z_dist_num(:,1).^2 + z_dist_num(:,2).^2 + z_dist_num(:,3).^2);
    
    z_dist = z_dist_num./z_dist_denom;
    
    z_dist = z_dist * side;
    
    %x-axis
    
    RFH_RLMA = [(PTC_x(:,b) - RFH_x(:,b)) (PTC_y(:,b) - RFH_y(:,b)) (PTC_z(:,b) - RFH_z(:,b))];
    
    x_dist_num = cross(RFH_RLMA,z_dist);

    %     RFH_RLMA = [(RFH_x(:,b) - RLMA_x(:,b)) (RFH_y(:,b) - RLMA_y(:,b)) (RFH_z(:,b) - RLMA_z(:,b))];
%     
%     x_dist_num = cross(RFH_RLMA,z_dist);

    
    x_dist_denom = sqrt(x_dist_num(:,1).^2 + x_dist_num(:,2).^2 + x_dist_num(:,3).^2);
    
    x_dist = x_dist_num./x_dist_denom;
    
    %y-axis
    
    y_dist = cross(z_dist,x_dist);

    %Setting up the Joint Coordinate System
    
    FA = cross(y_dist,z_prox); %Calculates the floating axis
    
    alpha(:,b) = 180 - acosd(dot(FA,y_prox,2)); %Holds the flexion 
    gamma(:,b) = asind(dot(FA,z_dist,2)); %Holds internal/external rotation
    beta(:,b) = acosd(dot(z_prox,y_dist,2)); %Holds abduction angle
     
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
    
    z_dist_2_num = [(RMP5_x(:,b) - R2M_x(:,b)) (RMP5_y(:,b) - R2M_y(:,b)) (RMP5_z(:,b) - R2M_z(:,b))];
    
    z_dist_2_denom = sqrt(z_dist_2_num(:,1).^2 + z_dist_2_num(:,2).^2 + z_dist_2_num(:,3).^2);
    
    z_dist_2 = z_dist_2_num./z_dist_2_denom;
    
    z_dist_2 = z_dist_2 * side;
    
    %x-axis
    
    RCAL_R5M = [(RCAL_x(:,b) - RMP5_x(:,b)) (RCAL_y(:,b) - RMP5_y(:,b)) (RCAL_z(:,b) - RMP5_z(:,b))];
    
    x_dist_2_num = cross(RCAL_R5M,z_dist_2);
    
    x_dist_2_denom = sqrt(x_dist_2_num(:,1).^2 + x_dist_2_num(:,2).^2 + x_dist_2_num(:,3).^2);
    
    x_dist_2 = x_dist_2_num./x_dist_2_denom;
    
    %y-axis
    
    y_dist_2 = cross(z_dist_2,x_dist_2);
    
    %Floating Axis
    
    FA_2 = cross(y_dist_2,z_prox_2);
    
    %Range of Motion of Joint
    
    alpha_2(:,b) = 180 - asind(dot(FA_2,y_prox_2,2)); %Holds the flexion 
    gamma_2(:,b) = asind(dot(FA_2,z_dist_2,2)); %Holds internal/external rotation
    beta_2(:,b) = acosd(dot(z_prox_2,y_dist_2,2)); %Holds abduction angle
    

% Not all segments have the same length for some reason which is causing
% prox_x from the BA to have an irregular length
%Establish Joint Coordinate System for BRACHIUM and SCAPULA
    
%Local Coordinate System for PROXIMAL segment
  
%z-axis

z_prox_3_num = [(RIWG_x(:,b) - RISC_x(:,b)) (RIWG_y(:,b) - RISC_y(:,b)) (RIWG_z(:,b) - RISC_z(:,b))];
    
z_prox_3_denom = sqrt(z_prox_3_num(:,1).^2 + z_prox_3_num(:,2).^2 + z_prox_3_num(:,3).^2);
    
z_prox_3 = z_prox_3_num./z_prox_3_denom;

z_prox_3 = z_prox_3 * side;
    
%x-axis
    
CRS_RISC = [(CRS_x(:,b) - RISC_x(:,b)) (CRS_y(:,b) - RISC_y(:,b)) (CRS_z(:,b) - RISC_z(:,b))];
    
    x_prox_3_num = cross(CRS_RISC,z_prox_3);
    
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
    
    %Range of Motion of joint (HIP)
      
    alpha_3(:,b) = 180 - asind(dot(FA_3,y_prox_3,2)); %Holds the flexion 
    gamma_3(:,b) = asind(dot(FA_3,z_dist_3,2)); %Holds internal/external rotation
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

