function [gc_is_good] = error_check(M1,M2,frame_errors,error_percent,static_seg,start_frame,end_frame)
    
    M1 = M1(start_frame:end_frame,:);
    M2 = M2(start_frame:end_frame,:);
    
    [M1_x, M1_y, M1_z] = extract_XYZ(M1);
    [M2_x, M2_y, M2_z] = extract_XYZ(M2);

    if(isnan(frame_errors))
        frame_errors = 15;
    end
    
%     if(isnan(error_percent))
%         error_percent = 0.03;
%     end
    
        
    error_count = 0;
    is_good = 1;

    %{
    %Constructs and compares anatomical segments for the static and dynamic
    %trials. The skeleton is assumed to be a rigid body, where the relative distance between
    %two points on the rigid body are always in constant distance to each
    %other: but because of skin motion and the potential for readings to be
    %missed by the camera's capture rate, the dynamic trial will produce results
    %where two given anatomical points are not within constant distance of each
    %other at all times. 

    %This section determines if the difference of the distance between two
    %points during the dynamic trial is too different from the static trial to
    %be able to conclude that the dynamic trial recorded accurately represents
    %the motion of the dog's skeleton.
    %}
   
       x = M1_x(1) - M2_x(1); 
       z = M1_z(1) - M2_z(1);
       y = sqrt(x^2 + z^2);
       
       k = 0; %current error count
       
       kk = 0;
       k_good = 1;
       
    for kk = 1:length(M1_x)

        AB_x = M1_x(kk) - M2_x(kk); 
        AB_z = M1_z(kk) - M2_z(kk);
        AB = sqrt(AB_x^2 + AB_z^2);
        
        low_accept = static_seg - error_percent;
        hi_accept = static_seg + error_percent;
        
%        if(mod(kk,2) == 0) 
%             disp("Current frame # " + kk);
%             disp("Dynamic: " + AB);
%             disp("Static: " + static_seg);
%             disp("Lower bound: " + low_accept);
%             disp("Upper bound: " + hi_accept);
%             disp("current count: "  + k);
%             disp("Actual difference between stat and dyn segments " + abs(AB - static_seg))
%        end
       
              
        if((AB <= low_accept) || (AB >= hi_accept))
            k = k + 1;
%             disp("Frame marked as error");
            error_count = error_count + 1;   
        else
            k = 0;
        end
        
        if(k >= frame_errors)
        
            k_good = 0;      
            break;
        end
        
     end %for
    
   if(k_good)
       gc_is_good = 1;
   elseif(~k_good)
       gc_is_good = 0;
   end
   
end %error_check




