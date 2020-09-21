function [alpha,gamma,beta,ps1,ps2,ps3] = create_angle_data(marker1,marker2,marker3,marker4,marker5,marker6,marker7,ps1,ps2,ps3)

   
    [marker1_x marker1_y marker1_z] = extract_XYZ(marker1);
    [marker2_x marker2_y marker2_z] = extract_XYZ(marker2);
    [marker3_x marker3_y marker3_z] = extract_XYZ(marker3);
    [marker4_x marker4_y marker4_z] = extract_XYZ(marker4);
    [marker5_x marker5_y marker5_z] = extract_XYZ(marker5);
    [marker6_x marker6_y marker6_z] = extract_XYZ(marker6);
    [marker7_x marker7_y marker7_z] = extract_XYZ(marker7);
    
    b = 1;
    
    if(marker1_x(1) ~= -1)
        %z-axis
        z_prox_num = [(marker1_x(:,b) - marker2_x(:,b)) (marker1_y(:,b) - marker2_y(:,b)) (marker1_z(:,b) - marker2_z(:,b))]; 
        z_prox_denom = sqrt(z_prox_num(:,1).^2 + z_prox_num(:,2).^2 + z_prox_num(:,3).^2);
        z_prox = z_prox_num./z_prox_denom;

        %x-axis 
        seg1 = [(marker3_x(:,b) - marker1_x(:,b)) (marker3_y(:,b) - marker1_y(:,b)) (marker3_z(:,b) - marker1_z(:,b))];   
        x_prox_num = cross(seg1,z_prox);
        x_prox_denom = sqrt(x_prox_num(:,1).^2 + x_prox_num(:,2).^2 + x_prox_num(:,3).^2);    
        x_prox = x_prox_num./x_prox_denom;

        %y-axis    
        y_prox = cross(z_prox,x_prox);
           
        else %handles AP calculations
            z_prox = ps1;
            x_prox = ps2;
            y_prox = ps3;
     end
        

    %Local Coordinate System for DISTAL segment      
    if(marker4_x(1) ~= -1 && marker5_x(1) ~= -1)
            
        %z-axis   
        z_dist_num = [(marker4_x(:,b) - marker5_x(:,b)) (marker4_y(:,b) - marker5_y(:,b)) (marker4_z(:,b) - marker5_z(:,b))];   
        z_dist_denom = sqrt(z_dist_num(:,1).^2 + z_dist_num(:,2).^2 + z_dist_num(:,3).^2);   
        z_dist = z_dist_num./z_dist_denom;

        %x-axis   
        seg2 = [(marker6_x(:,b) - marker4_x(:,b)) (marker6_y(:,b) - marker4_y(:,b)) (marker6_z(:,b) - marker4_z(:,b))];   

        x_dist_num = cross(seg2,z_dist);   
        x_dist_denom = sqrt(x_dist_num(:,1).^2 + x_dist_num(:,2).^2 + x_dist_num(:,3).^2);
        x_dist = x_dist_num./x_dist_denom;

        %y-axis  
        y_dist = cross(z_dist,x_dist);

        %Setting up the Joint Coordinate System
        FA = cross(y_dist,z_prox); %Calculates the floating axis
    else %BS case
        z_dist = ps1;
        x_dist = ps2;
        y_dist = ps3;
        FA = cross(z_prox,y_dist);
    end

        alpha(:,b) = 180 - acosd(dot(FA,x_prox,2)); %Holds the flexion 
        gamma(:,b) = acosd(dot(z_dist,FA,2)); %Holds internal/external rotation
        beta(:,b) = acosd(dot(z_prox,y_dist,2)); %Holds abduction angle
        
        if(marker1_x ~= -1)
            ps1 = z_dist;
            ps2 = x_dist;
            ps3 = y_dist;
        else %AP case
            ps1 = z_prox;
            ps2 = x_prox;
            ps3 = y_prox;
        end
end
