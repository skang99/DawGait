function [gait_cycle_loc, gait_cycle_count] = find_gait_cycle_frames(R_5th_M_z)


test_tol = mean(R_5th_M_z);

gait_cycle_loc = 1;

%This section finds the locations of where each gait cycle begins (aside
%from the first gait cycle, whose location was calculated previously)

abs_max = max(R_5th_M_z); %finds the absolute maximum position of the new Z-paw array

abs_min = min(R_5th_M_z); %finds the absolute minimum position of the new Z_paw array

max_min_tol = (abs_min + abs_max)/2; %sets a tolerance to help find the beginning of each gait cycle

new_max_data = findpeaks(R_5th_M_z); %finds all local maxima in the new Z-paw array
                                     %and stores it in a seperate array

new_min_data = findpeaks(-R_5th_M_z)*-1; %finds all local minima in the Z-paw array
                                         %and stores it in a seperate array
                                                                      
count = 0; 

for xx = 1:length(new_min_data)     
    max_min_diff = new_max_data(xx) -  new_min_data(xx);
    
    if max_min_diff >= test_tol        
        count = count + 1; %Keeps track of how many gait cycles have been found                          
                                   
        new_min_pos(count) = new_min_data(xx); %Finds position where the gait cycle begins
          
    end
    
end %ends loop

gait_cycle_count = count;


%count is equal to number of gait cycles - 1
    if count == 0

        gait_cycle_loc = -1;

    else

        gait_cycle_loc(1) = 1;

        for i = 1:length(new_min_pos)

            t = find(R_5th_M_z == new_min_pos(i));

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

         gait_cycle_loc(i+1) = t;
            
        end  %for

    end
    

    
end