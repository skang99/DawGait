%gait cycle averaging function
clearvars;
clc;
add = true;
tobeavgd = [];
jointnames = ["SHLD" "ELB" "CARP"];
dirnames = ["_f" "_r" "_a"];
%select joint
joint = input("Which Joint? 1-Shoulder 2-Elbow 3-Wrist \n");
%select dimension
dir = input("Which dimension? 1-Flexion 2-Rotation 3-Abduction \n");
specifier = strcat(jointnames(joint),dirnames(dir));
%load file
while add
    [file,path] = uigetfile();
    if file ~= 0
        load(strcat(path,file));
    else
        disp("Program execution cancelled")
        return
    end
    %normalize to 100
    xx = 1:100;
    j = fieldnames(jc_angles.angles)';
    j = j(1:9);
    for i = 1:size(jc_angles.angles,2)
        for k = 1:9
            tseries = jc_angles.angles(i).(j{k});
            if tseries(1) ~= -1
                yy = spline(1:length(tseries),tseries,linspace(1,length(tseries)));
                jc_angles.angles(i).(j{k}) = yy';
                %disp("done");
            else
                jc_angles.angles(i).(j{k}) = (-1 * ones(100,1));
                %disp("bad gc")
            end
        end
    end
    %select gait cycles
    cycles = input(sprintf("enter gait cycle number(s) to be included in square brackets separated by spaces(trial has %d gait cycles): \n",size(jc_angles.angles,2)));
    %add them to a list of gait cycles to average
    for i = 1:length(cycles)
        tobeavgd = [tobeavgd jc_angles.angles(cycles(i)).(specifier)];
    end
    %ask if more files to add
    cont = input("Add more files? y/n: \n",'s');
    if cont == 'y' || cont == 'Y'
        add = true;
    else
        add = false;
    end
end
%compute the averages
tobeavgd
numfiles = size(tobeavgd,2);
GaitCycleAvg = sum(tobeavgd,2)
GaitCycleAvg = GaitCycleAvg./numfiles;
%export to file
save('Gait Cycle Averages','GaitCycleAvg','tobeavgd');



