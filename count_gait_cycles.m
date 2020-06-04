function [gait_cycle_loc] = count_gait_cycles(R_5th_M_z)

%This function finds the specific starting and ending frames of each gait cycle.
%Unlike find_cycle_frames, this function will return an array of the
%specific starts and ends of each gait cycle instead of just the start of
%the first and end of the last cycles. The precondition for this function
%is that the position data has been trimmed to only include the position
%data found for the gait cycle.

%Sets a tolerance to help find the beginning of each gait cycle.
gc_tolerance = mean(R_5th_M_z);

new_max_data = findpeaks(R_5th_M_z);
new_min_data = findpeaks(-R_5th_M_z)*-1;

%The start of a gait cycle is found by determining if a given local
%minima's relative distance to a corresponding local maxima is greater than
%the average of all gait cycle position data. Because gait cycles are
%determined to begin and end before and after some apex point, this ensures
%that the start of a gait cycle is the dog's foot reaching the ground and
%then lifting fully.

count = 0; 

for xx = 1:length(new_min_data) 
    max_min_diff = new_max_data(xx) -  new_min_data(xx) 
    
    if max_min_diff >= gc_tolerance %deterines whether or not the difference previously calculated
                                   %is greater than or equal to average of the absolute maximum and absolute minimum
        
        count = count + 1; %Keeps track of how many gait cycles have been found                          
                                   
        new_min_pos(count) = new_min_data(xx); %Finds position where the gait cycle begins
          
    end %if
end %for

%Because the position data for the gait cycle has already been trimmed, if
%no other acceptable start points have been found, only one gait cycle has
%been found.
if count == 0
    gait_cycle_loc = -1;
    return;
else
    gait_cycle_loc(1) = 1;
    for i = 1:length(new_min_pos) 
        t = find(R_5th_M_z == new_min_pos(i))   
        gait_cycle_loc(i+1) = t;         
    end %for  
end %if

end

