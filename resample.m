function [] = resample(trial)

load(trial);

xx = 1:100;
jc_angles_resampled = jc_angles;
j = fieldnames(jc_angles.angles)';
j = j(1:9);
for i = 1:size(jc_angles.angles,2)
    for k = 1:9
        tseries = jc_angles.angles(i).(j{k});
        if tseries(1) ~= -1
            yy = spline(1:length(tseries),tseries,linspace(1,length(tseries)));
            jc_angles_resampled.angles(i).(j{k}) = yy';
            disp("done");
        else
            jc_angles_resampled.angles(i).(j{k}) = (-1 * ones(100,1));
            disp("bad gc")
        end
    end
end
newname = jc_angles.trial_name + " Angles Resampled.mat";
save(newname,'jc_angles_resampled');


end
