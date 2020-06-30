function [gait_cycle_loc, gait_cycle_count] = find_gait_cycle_frames(R_5th_M_z)

test_tol = mean(R_5th_M_z)

gait_cycle_loc = 1;

%This section finds the locations of where each gait cycle begins (aside
%from the first gait cycle, whose location was calculated previously)

abs_max = max(R_5th_M_z); 

abs_min = min(R_5th_M_z); 

max_min_tol = (abs_min + abs_max)/2; %sets a tolerance to help find the beginning of each gait cycle


new_max_data = findpeaks(R_5th_M_z)

new_min_data = findpeaks(-R_5th_M_z)*-1;
                                         
%Empty positional readings trimmed
new_min_data = new_min_data(new_min_data ~= 0)

%Trims mins/maxes from the respective vectors that are above/below
%the average position: these points are certainly not the trough/peak of a gc

new_min_data = new_min_data(new_min_data < test_tol)
new_max_data = new_max_data(new_max_data > test_tol)
                                                                                                                    
count = 0; 

length(new_min_data)
length(new_max_data)

min_array = [];
max_array = [];

for i=1:length(new_min_data)
    min_array(i) = find(R_5th_M_z == new_min_data(i));
end

min_array

for i=1:length(new_max_data)
    max_array(i) = find(R_5th_M_z == new_max_data(i));
end

max_array

% Must only run n amount of checks where n is the length of the smaller of
% the min/max arrays

check_length = 0;

if(length(new_min_data) > length(new_max_data))
    check_length = length(new_max_data);
else
    check_length = length(new_min_data);
end

disp("min: " + length(new_min_data))
disp("max " + length(new_max_data))
disp(check_length)
    

%Finds gait cycles by checking the distance between a local max and min. If
%this value is greater than the tolerance (the average of the position
%data), the two points are likely far enough apart to be a peak and a
%trough of a gait cycle
for xx = 1:check_length
    
    max_min_diff = new_max_data(xx) -  new_min_data(xx)
    
    t_min = find(R_5th_M_z == new_min_data(xx))
    t_max = find(R_5th_M_z == new_max_data(xx))
    
 
    if max_min_diff > test_tol
       
        count = count + 1; %Keeps track of how many gait cycles have been found
        
        disp("success min " + new_min_data(xx))
        disp("success max " + new_max_data(xx))
        
        new_min_pos(count) = new_min_data(xx); %Stores an array of local mins that are determined to be the trough of a gait cycle
        
     
        %Iteration 5:
        %Local min at pos 32.6061 is determined to be a suitable local min
        %for after the previous gait cycle's end at pos 24.6104
        %since this iteration's local max at 89.9418 falls within tolerance
        
        %Iteration 6:
        %local max at 95.2338 corresponds to the local min at 23.7633 so
        %that min is also counted as the end of a gait cycle
    end
    
end %ends loop

gait_cycle_count = count;



%count is equal to number of gait cycles - 1
    if count == 0

        gait_cycle_loc = -1;

    else
          
        gait_cycle_loc(1) = 1;

        for i = 1:length(new_min_pos)
            
            disp("checking " + new_min_pos(i))

            t = find(R_5th_M_z == new_min_pos(i)) %finds the specific frame the local min occurs

            %Fix for cases where the calculated starting position value within tolerance for
            %a gait cycle's beginning is found at multiple times. Assuming
            %the dog moves with a somewhat consistent speed through the
            %gait cycles, the gait cycles are roughly spaced equidistantly
            %through the trial duration.
            
            if(length(t) > 1)  
                tolerance_frame = length(R_5th_M_z) / (gait_cycle_count + 1);

                 for j = 1:length(t)
                     test_frame = t(j);                
                     frame_distances(j) = abs(test_frame - tolerance_frame);
                 end

                 best_tolerance = min(frame_distances);

                 for j = 1:length(frame_distances)
                     if(frame_distances(j) == best_tolerance) 
                         t = t(j);
                         break;
                     end %if
                 end %for

            end %if

            %Assumes the next local min is the end of a gait cycle and the
            %start of the next one
            if(i ~=1)
                disp("forward diff: " + abs(gait_cycle_loc(i-1) - t))
            end
            
        gait_cycle_loc(i+1) = t;
                
            
        gait_cycle_loc
            
        end  %for

    end %if
    
    

    
end

