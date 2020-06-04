function [t_naught,t_one] = find_still_static(R_5th_M_z)
% Finds the period of time in which the dog is the most still during
% the static trial.

frames_per_second = 200;
means = zeros(1,floor(length(R_5th_M_z) / frames_per_second));
second = 1;
sum = 0;

%Finds the average of the z-position for every one second, 100 frames, 
%and selects the smallest average to be the period in which the dog is
%still

for i = 1:length(R_5th_M_z)
    
    sum = sum + R_5th_M_z(i);
    
    if(mod(i,frames_per_second) == 0)
        means(second) = sum/frames_per_second;
        second = second + 1;
    end
    
end

min_average = min(means);

t_one = find(means == min_average);
t_naught = t_one - 1;
    
  
end

