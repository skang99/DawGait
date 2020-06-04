%returns the frames at which the gait cycle(s) begin and end
function [cyc_start,cyc_end] = find_cycle_frames(marker_pos)


% The purpose of this function is to find the frame at which the first gait
% cycle begins by checking if local minima of a given landmark's position graph
% are sufficiently low enough, ie, below an average position, to be
% considered the beginning of the gait cycle.

max_pos = max(marker_pos);
min_pos = min(marker_pos);
avg_pos = (max_pos + min_pos)/2;
min_data = findpeaks(-marker_pos)*-1; %finds all local minima in the z-paw data and puts it into an array

for j = 1:length(min_data) %starts loop to find the first suitable minimum
    
    if min_data(j) < avg_pos %determines whether the local minimum is less than the specified tolerance
        
        first_min_pos = min_data(j); %A z coordinate that will correspond to a frame number at which the beginning of the gait cycle starts
        
        break; %breaks from loop
    end % end if statement
end %ends loop

first_min_loc = find(marker_pos == first_min_pos); %finds the frame for the corresponding minimum value in the Z-paw array

cyc_start = first_min_loc;

%Once the exact frame is found, because all anatomical landmarks move at
%the same time, the gait cycle's beginning in each landmarks' position
%graph can be found.

%Everything before the start of the full gait cycle is deleted
marker_pos(1:first_min_loc-1) = [];

% This section finds the frame where the final, full gait cycle ends

rev_z_paw = flipud(marker_pos); %flips the z-paw data up to down
r2 = mean(rev_z_paw(1:10)); %calculates the average position of the last 10 z-paw positions

if r2 >= avg_pos
    max_tol = r2 + 0.05*r2; 
else  
    max_tol = r2 + 0.2*r2; 
end

max_data = findpeaks(marker_pos); %finds all local maxima within the Z-paw array
rev_max_data = flipud(max_data); %reverses the array containing all local maxima in order to find the last suitable max position of the z-paw

for z = 1:length(rev_max_data) %starts loop to find the max value we need
    if rev_max_data(z) >= max_tol %determines whether the local maximum is greater than or equal to the sepecified tolerance. If it is not then the loop continues. 
        max_pos = rev_max_data(z); %finds the maximum position of the Z-paw we need     
        break; %breaks from loop
    end %ends if statement
end %ends loop

max_loc = find(marker_pos == max_pos); %finds the frame at which the last suitable maximum value is located

second_min_data = findpeaks(-marker_pos(max_loc:end))*-1; %finds all local minimum that follow the max position
                                                         %calculated in the previous loop within the Z-paw array
                                                         
second_min_pos = second_min_data(1); %finds the position of the Z-paw where the final, complete gait cycle ends

second_min_loc = find(marker_pos == second_min_pos); %finds the frame where the final, complete gait cycle ends

marker_pos(second_min_loc:end,:) = [];

%The purpose of this section is remove data from the extracted position
%readings that do not form a complete gait cycle: only data that produces a
%full gait cycle will be kept.

a = max(marker_pos);
b = min(marker_pos);
c = a - b;
tol = c - c*0.2;

max_data = findpeaks(marker_pos);
d = max_data(end) - marker_pos(end);

if d < tol
       
    min_data = findpeaks(-marker_pos)*-1;
    
    max_min_tol = (a + b)/2;
    
    max_min_tol_2 = max_min_tol - max_min_tol*0.1;
    
    for k = 1:length(max_data)
      
        max_min_diff = max_data(k) - min_data(k);
        
        if max_min_diff >= max_min_tol_2
            
            min_pos = min_data(k);
            
            break;
        end
        
    end   
    
    min_loc = find(marker_pos == min_pos);
    marker_pos(min_loc:end,:) = [];
   
end %if

cyc_end = cyc_start+length(marker_pos);

end

