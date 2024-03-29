function [cyc_start,cyc_end] = find_start_cycle_frame(marker_pos)

% The purpose of this function is to find the frame at which the first gait
% cycle begins and the last gait cycle ends, defined by local minima
% marker_pos may not be 0 or blank

max_pos = max(marker_pos);
min_pos = min(marker_pos);
avg_pos = (max_pos + min_pos)/2;
min_data = findpeaks(-marker_pos)*-1;
max_data = findpeaks(marker_pos);

min_data = min_data(min_data < (avg_pos - avg_pos * 0.15) & min_data ~= 0);
max_data = max_data(max_data > avg_pos);

first_min_pos = min_data(1);

for j = 1:length(min_data) 
    %In order for a min to be considered the start of a gait cycle, the min
    %must be less than the average of the position data and is defined as
    %the point where the dog's foot is on the ground
    
    if (min_data(j) < avg_pos)  
        first_min_pos = min_data(j); %A z coordinate that will correspond to a frame number at which the beginning of the gait cycle start 
        break;         
    end % if
    
end %for

first_min_loc = find(marker_pos == first_min_pos,1); %finds the frame for the corresponding minimum value in the Z-paw array

cyc_start = first_min_loc;

%Once the exact frame is found, because all anatomical landmarks move at
%the same time, the gait cycle's beginning in each landmarks' position
%graph can be found.

%Everything before the start of the full gait cycle is deleted
marker_pos(1:first_min_loc-1) = [];

% This section finds the frame where the final, full gait cycle ends
rev_min_data = flipud(min_data);
final_min = rev_min_data(1);
    
min_loc = find(marker_pos == final_min,1);
marker_pos(min_loc:end,:) = [];

cyc_end = cyc_start+length(marker_pos);
    
   
end

