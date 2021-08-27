function [resampled_angles,gc_count,gc_frames,trial_angles] = resample_trial(trial,joint,plane)

joint_names = ["SHLD" "ELB" "CARP"];
plane_names = ["_f" "_r" "_a"];

specifier = strcat(joint_names(joint),plane_names(plane));

load(['Produced Data/' trial]);
total_gc_count = length(jc_angles.angles);
% Counts GCs with non blank data
gc_count = 0;

k = 1;
for i = 1:total_gc_count
    
    % Removes the frames for GCs with blank data
    if(jc_angles.angles(i).(specifier)(1) ~= -1)
        gc_count = gc_count + 1;
        gc_frames(k) = jc_angles.angles(i).start_frame;
        gc_frames(k+1) = jc_angles.angles(i).end_frame;
        k = k+1;
    end
end

%Readjusts GC locations to begin at frame #1 and end at cyc_end
if(gc_frames(1) ~= 1)  
    gc_frames = gc_frames - gc_frames(1);
    gc_frames(1) = 1;
end

trial_angles = jc_angles.(specifier);

j = fieldnames(jc_angles.angles)';
% Truncates the start_frame and end_frame field names
% {'ELB_f','ELB_r','ELB_a','CARP_f','CARP_r','CARP_a','SHLD_f','SHLD_r','SHLD_a'}
j = j(1:9);

resampled_angles = [];

for i = 1:total_gc_count
    tseries = jc_angles.angles(i).(specifier);
    if tseries(1) ~= -1
        yy = spline(1:length(tseries),tseries,linspace(1,length(tseries)));
        resampled_angles = [resampled_angles; yy];
    end
end

% Ensures trial data length agrees with dimension of trial's duration
trial_angles = trial_angles(1:gc_frames(length(gc_frames)));

resampled_angles = (sum(resampled_angles',2)./gc_count)';


end
