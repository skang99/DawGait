function [cycle_frame, gait_cycle_count] = split_gait_cycles(marker_pos)

% Assumes marker_pos has been trimmed of leading and trailing data: this
% function executes after find_start_cycle_frame finds the start of the
% first gait cycle and end of the last gait cycle: therefore, the frame
% numbers returned by this function will not line up with the frame numbers of
% the recorded z positional data in the dynamic trial. For clarity, the gait
% cycles frames will adjusted outside of this function to begin at 1 before
% being placed in the mat file, where that frame will represent the start of 
% all gait cycles and the start of gait cycle 1

% Note that marker_pos should be z data. Trials with poorly recorded z data
% for the 5th metacarpal (our standard landmark for a trial's recording quality) will produce 
% incorrect results and should be discarded: blank data or spotty data,
% particularly near frames approaching bowls, will make this approach unreliable

max_data = findpeaks(marker_pos);
min_data = findpeaks(-marker_pos) * -1;
avg_m = mean(min_data);

%trimmed pos data
avg = mean(marker_pos);
abs_min = min(marker_pos);

%Trims maxes/mins that are not above/below the mean of data
max_data = max_data(max_data > avg);
min_data = min_data(min_data < avg & min_data ~= 0);

%Trims maxes/mins that are not above/below the mean of data
max_data = max_data(max_data > avg);
min_data = min_data(min_data < avg & min_data ~= 0);

max_frames = [];
min_frames = [];

% Records frame numbers for "filtered" maxima that are actually above avg
for i=1:length(max_data)
    max_frames(i) = find(marker_pos == max_data(i),1);
end

%absolute min - mean of min_data
e = abs((abs_min - avg_m)) / avg_m;

% Tolerance value
avg_e = avg - (avg * e);

min_data = min_data(min_data < avg_e);

% Records frame numbers for "filtered" minima that are actually below avg
for i=1:length(min_data)
    min_frames(i) = find(marker_pos == min_data(i),1);
end

p_sign = marker_pos(1) > 0;
n_sign = marker_pos(2) > 0;
index = 1;
bowl_test = [];

% Bowl testing: the area between 2 local maxima makes a bowl where avg is
% the top of the bowl

% Finds intersects, "bowl" starts and ends
% Determines the frame where marker_pos changes from being below to above
% the average of all data
for i=1:length(marker_pos)
    if(marker_pos(i) - avg > 0)
        n_sign = 1;
    elseif(marker_pos(i) - avg < 0)
        n_sign = 0;
    else
        bowl_test(index) = i;
        index = index+1;
    end %if
    
    if(p_sign ~= n_sign)
        bowl_test(index) = i;
        index = index+1;
    end
    
    % Compared again in the next iteration
    p_sign = n_sign;
end

%The first bowl should start at 1 and the end of the final bowl should end
%at the last pos value
bowl_test(length(bowl_test) + 1) = length(marker_pos);

index = 1;
check_frames = [];
pos_data = [];

% Finds gait cycle frames
for i=1:length(bowl_test) - 1
    z1 = bowl_test(i);
    z2 = bowl_test(i+1);
        
    %Tests bowl frames two at a time, checks if the two tested frames
    %contain local minima found before, as points that intersected the avg
    %value are not necessarily bowls unless they contain the local mins
    %that were found to be below the avg of marker_pos
    for k = 1:length(min_frames)
        if(min_frames(k) > z1 & min_frames(k) < z2)
            check_frames(k) = min_frames(k); %holds frames of mins in this interval (z1,z2) to be checked
            pos_data(k) = marker_pos(min_frames(k)); %holds pos data of those frames
        end
    end %for 
    
    break_loop = 0;

    % A bowl should not have a filtered maxima inside of it: if this
    % condition is not met, the frames z1 and z2 do not make a correct bowl
    % between two gait cycles 
    for h = 1:length(max_frames)
        if(max_frames(h) > z1 & max_frames(h) < z2)
            break_loop = 1;
            break
        end
    end %for
    
    if(break_loop)
        continue
    end %if
    
    %z1 and z2 are now defined a bowl
    bowl_pos_data = marker_pos(z1:z2);

    current_min = bowl_pos_data(1);
    
    %The smallest value in the bowl will then be the end of the gait cycle
    %and the start of the next 
    for q=1:length(bowl_pos_data)
        if(current_min > bowl_pos_data(q))
            current_min = bowl_pos_data(q);
        end
    end %for
    
    lowest_frame = find(marker_pos == current_min,1);
    cycle_frame(index) = lowest_frame;
    index = index + 1;
                 
    %Trims the min_frame array to not include frames of mins already
    %checked
    min_frames = setdiff(min_frames,check_frames);
    check_frames = [];      
    pos_data = [];
   
end

gait_cycle_count = length(cycle_frame) - 1;

if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
end

end





    
    
    
            






