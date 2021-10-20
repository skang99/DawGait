function [gc_is_good] = error_check(M1,M2,error_length,static_seg,start_frame,end_frame,t)
% Returns whether the segment created by M1 and M2 passes the error check
% One or both markers may be empty: all blank data was understood to be 0

M1 = M1(start_frame:end_frame,:);
M2 = M2(start_frame:end_frame,:);

[M1_x, M1_y, M1_z] = extract_XYZ(M1);
[M2_x, M2_y, M2_z] = extract_XYZ(M2);

% Acceptable number of consecutive error frames
frame_errors = 15;

if(isnan(error_length))
    error_length = 50; %mm
end

error_count = 0;

%{
    %Constructs and compares anatomical segments for the static and dynamic
    %trials. The skeleton is assumed to be a rigid body, where the relative distance between
    %two points on the rigid body are always in constant distance to each
    %other: because of skin motion and the potential for readings to be
    %missed by the camera's capture rate, the dynamic trial will produce results
    %where two given anatomical points are not within constant distance of each
    %other at all times.

    %This section determines if the difference of the distance between two
    %points during the dynamic trial is too different from the static trial to
    %be able to conclude that the dynamic trial recorded accurately represents
    %the motion of the dog's skeleton.
%}

k = 0; %current error count

k_good = 1;

for kk = 1:length(M1_x)
    
    AB_x = M1_x(kk) - M2_x(kk);
    AB_z = M1_z(kk) - M2_z(kk);
    AB = sqrt(AB_x^2 + AB_z^2);
    
    low_accept = static_seg - error_length;
    hi_accept = static_seg + error_length;
    
    % low_accept <= AB <= hi_accept
    if((AB <= low_accept) || (AB >= hi_accept))
        k = k + 1;
        error_count = error_count + 1;
    else
        k = 0;
    end %if
    
    if(k >= frame_errors)
        k_good = 0;
        break;
    end %if
    
end %for

if(k_good)
    gc_is_good = 1;
elseif(~k_good)
    gc_is_good = 0;
end


end %error_check




