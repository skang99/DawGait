function [cycle_frame, gait_cycle_count] = split_gait_cycles(marker_pos)

% avg = mean(marker_pos);
% abs_min = min(marker_pos);
 

max_data = findpeaks(marker_pos);
min_data = findpeaks(-marker_pos) * -1;
avg_m = mean(min_data);

%trimmed pos data
avg = mean(marker_pos);
abs_min = min(marker_pos);

%Trims maxes/mins that are not above/below the mean of data
max_data = max_data(max_data > avg);
min_data = min_data(min_data < avg & min_data ~= 0);

% %trimmed pos data
% avg = mean(marker_pos);
% abs_min = min(marker_pos);

%Trims maxes/mins that are not above/below the mean of data
max_data = max_data(max_data > avg);
min_data = min_data(min_data < avg & min_data ~= 0);

max_frames = [];
min_frames = [];

for i=1:length(max_data)
    max_frames(i) = find(marker_pos == max_data(i),1);
end

e = abs((abs_min - avg_m)) / avg_m;
avg_e = avg - (avg * e);

min_data = min_data(min_data < avg_e);

for i=1:length(min_data)
    min_frames(i) = find(marker_pos == min_data(i),1);
end

p_sign = marker_pos(1) > 0;
n_sign = marker_pos(2) > 0;
index = 1;
bowl_test = [];

%Finds intersects, "bowl" starts and ends
for i=1:length(marker_pos)
    if(marker_pos(i) - avg > 0)
        n_sign = 1;
    elseif(marker_pos(i) - avg < 0)
        n_sign = 0;
    else
        bowl_test(index) = i;
        index = index+1;
    end
    
    if(p_sign ~= n_sign)
        bowl_test(index) = i;
        index = index+1;
    end
    
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
    %contain minimums found before, as points that intersected the avg
    %value are not necessarily actual bowls unless they contain mins
    for k = 1:length(min_frames)
        if(min_frames(k) > z1 & min_frames(k) < z2)
            check_frames(k) = min_frames(k); %holds frames of mins in this interval (z1,z2) to be checked
            pos_data(k) = marker_pos(min_frames(k)); %holds pos data of those frames
        end
    end
    
    break_loop = 0;

    for h = 1:length(max_frames)
        if(max_frames(h) > z1 & max_frames(h) < z2)
            break_loop = 1;
            break
        end
    end
    
    if(break_loop)
        continue
    end

    %check all pos data approach
    bowl_pos_data = marker_pos(z1:z2);

    current_min = bowl_pos_data(1);
    
    for q=1:length(bowl_pos_data)
        if(current_min > bowl_pos_data(q))
            current_min = bowl_pos_data(q);
        end
    end
    
    lowest_frame = find(marker_pos == current_min,1);
    cycle_frame(index) = lowest_frame;
    index = index + 1;
                 
    %Trims the min_frame array to not include frames of mins already
    %checked
    min_frames = setdiff(min_frames,check_frames);
    check_frames = [];      
    pos_data = [];
    d = [];
    
end

gait_cycle_count = length(cycle_frame) - 1;

if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
end

end





    
    
    
            






