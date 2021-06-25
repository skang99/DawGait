function [BA_f,BA_r,BA_a,AP_f,AP_r,AP_a,BS_f,BS_r,BS_a,new_time,R_5th_M_z,frames,no_good_cycles,gait_cycle_count] = hindlimbJCS(dynamic_trial,static_trial)

load(['Produced Data/' dynamic_trial])

gait_cycle_count = length(trial.gait_cycles);

no_good_cycles = 0;

MP5 = trial.pos_data.MP5;
[M5_x M5_y M5_z] = extract_XYZ(MP5);

MP2 = trial.pos_data.MP2;
LEP = trial.pos_data.LEP;
QUA = trial.pos_data.QUA;
MEP = trial.pos_data.MEP;
MMA = trial.pos_data.MMA;
GT = trial.pos_data.GT;
FH = trial.pos_data.FH;
LMA = trial.pos_data.LMA;
CRS = trial.pos_data.CRS;
IWG = trial.pos_data.IWG;
GAS = trial.pos_data.GAS;
CAL = trial.pos_data.CAL;
ISC = trial.pos_data.ISC;

direction = 1;

if(strcmp(trial.direction,"left"))
    direction = -1;
end


end_frame = length(trial.pos_data.MP5);

new_time = 1:end_frame;

[BA_f BA_r BA_a BS_f BS_r BS_a AP_f AP_r AP_a] = create_hindlimb_angle_data(LEP,MEP,GT,LMA,MMA,FH,MP5,MP2,CAL,ISC,CRS,IWG,direction);


%Graphs the velocity of T1 during gc1

% s = trial.gait_cycles(1).start_frame;
% e = trial.gait_cycles(1).end_frame;

% vT1 = T1_y(s:e);
% vel = [];
% 
% 
% for i=1:length(vT1)-1
%     vel(i) = abs(((vT1(i+1) - vT1(i)) / 0.02) / 1000);
% end

% figure(9)
% plot((s:e-1)/100,vel)
% xlabel('seconds(s)')
% ylabel('velocity(m/s)')
% title('T1_y velocity during gc#1')

%Show anyway? 
if(no_good_cycles)
    BA_f =0;BA_r=0;BA_a=0;AP_f=0;AP_r=0;AP_a=0;BS_f=0;BS_r=0;BS_a=0;new_time=0; R_5th_M_z=0;frames=0;
    no_good_cycles = 1;
    return;
end

nBA_f = [];
nBA_r = [];
nBA_a = [];

nAP_f = [];
nAP_r = [];
nAP_a = [];

nBS_f = [];
nBS_r = [];
nBS_a = [];

for i = 1:gait_cycle_count
    
%     if(~trial.gait_cycles(i).is_good)
%         continue;
%     end
    
    
    start_frame = trial.gait_cycles(i).start_frame;
    end_frame = trial.gait_cycles(i).end_frame;
    gc_length = abs(start_frame - end_frame) + 1;
    
    if(trial.gait_cycles(i).GT_LEP && trial.gait_cycles(i).FH_LMA)
        disp("Creating knee data for gc # " + i)
        tBA_f = BA_f(start_frame:end_frame);
        tBA_r = BA_r(start_frame:end_frame);
        tBA_a = BA_a(start_frame:end_frame);
        nBA_f = [nBA_f; tBA_f];
        nBA_r = [nBA_r; tBA_r];
        nBA_a = [nBA_a; tBA_a];
    else
        tBA_f = zeros(gc_length,1) + -1;
        tBA_r = zeros(gc_length,1) + -1;
        tBA_a = zeros(gc_length,1) + -1;
    end
    
    if(trial.gait_cycles(i).CAL_MP5 && trial.gait_cycles(i).FH_LMA)
        disp("Creating ankle data for gc # " + i)
        tAP_f = AP_f(start_frame:end_frame);
        tAP_r = AP_r(start_frame:end_frame);
        tAP_a = AP_a(start_frame:end_frame);
        nAP_f = [nAP_f; tAP_f];
        nAP_r = [nAP_r; tAP_r];
        nAP_a = [nAP_a; tAP_a];
    else
        tAP_f = zeros(gc_length,1) + -1;
        tAP_r = zeros(gc_length,1) + -1;
        tAP_a = zeros(gc_length,1) + -1;
    end 
    
    if(trial.gait_cycles(i).GT_LEP && trial.gait_cycles(i).IWG_ISC && trial.gait_cycles(i).CRS_ISC)
        disp("Creating hip data for gc # " + i)
        tBS_f = BS_f(start_frame:end_frame); 
        tBS_r = BS_r(start_frame:end_frame);
        tBS_a = BS_a(start_frame:end_frame);
        nBS_f = [nBS_f; tBS_f];
        nBS_r = [nBS_r; tBS_r];
        nBS_a = [nBS_a; tBS_a];
    else
        tBS_f = zeros(gc_length,1) + -1;
        tBS_r = zeros(gc_length,1) + -1;
        tBS_a = zeros(gc_length,1) + -1;
    end
       
    %postcondition: all data included, excluding bad(seg check failure)/missing data
    x(i) = struct("KNEE_f",tBA_f,"KNEE_r",tBA_r,"KNEE_a",tBA_a,"ANKLE_f",tAP_f,"ANKLE_r",tAP_r,"ANKLE_a",tAP_a,"HIP_f",tBS_f,"HIP_r",tBS_r,"HIP_a",tBS_a,"start_frame",start_frame,"end_frame",end_frame);

end
frames = [];

for i = 1:length(trial.gait_cycles)
    frames(i) = trial.gait_cycles(i).start_frame;
    frames(i+1) = trial.gait_cycles(i).end_frame;
end

jc_angles = struct("trial_name",trial.trial_name,"angles",x,"KNEE_f",nBA_f,"KNEE_r",nBA_r,"KNEE_a",nBA_a,"ANKLE_f",nAP_f,"ANKLE_r",nAP_r,"ANKLE_a",nAP_a,"HIP_f",nBS_f,"HIP_r",nBS_r,"HIP_a",nBS_a);

temp_name = char(trial.trial_name + " Angles");

save(['Produced Data/' temp_name '.mat'],'jc_angles')

R_5th_M_z = M5_z;

name = {trial.trial_name
        ''
        'GC#'};

header = {'1' '' '' '' '' '' '' '' ''};
label = {'KNEE_f' 'KNEE_r' 'KNEE_a' 'ANKLE_f' 'ANKLE_r' 'ANKLE_a' 'HIP_f' 'HIP_r' 'HIP_a'};

filename = "Produced Data/" + trial.trial_name + " Angular Data.xlsx";

header = repmat(header,1,gait_cycle_count);
label = repmat(label,1,gait_cycle_count);
index = 10;

for i=1:gait_cycle_count - 1
    header(index) = num2cell(i + 1);
    index = index + 9;
end

writecell(name,filename,'Range','A1')
writecell(header,filename,'Range','B3');
writecell(label,filename,'Range','B4');

G1 = table(jc_angles.angles(1).KNEE_f,jc_angles.angles(1).KNEE_r,jc_angles.angles(1).KNEE_a,jc_angles.angles(1).ANKLE_f,jc_angles.angles(1).ANKLE_r,jc_angles.angles(1).ANKLE_a,jc_angles.angles(1).HIP_f,jc_angles.angles(1).HIP_r,jc_angles.angles(1).HIP_a);

if(gait_cycle_count > 1)
    G2 = table(jc_angles.angles(2).KNEE_f,jc_angles.angles(2).KNEE_r,jc_angles.angles(2).KNEE_a,jc_angles.angles(2).ANKLE_f,jc_angles.angles(2).ANKLE_r,jc_angles.angles(2).ANKLE_a,jc_angles.angles(2).HIP_f,jc_angles.angles(2).HIP_r,jc_angles.angles(2).HIP_a);
end

if(gait_cycle_count > 2)
   G3 = table(jc_angles.angles(3).KNEE_f,jc_angles.angles(3).KNEE_r,jc_angles.angles(3).KNEE_a,jc_angles.angles(3).ANKLE_f,jc_angles.angles(3).ANKLE_r,jc_angles.angles(3).ANKLE_a,jc_angles.angles(3).HIP_f,jc_angles.angles(3).HIP_r,jc_angles.angles(3).HIP_a);
end

if(gait_cycle_count > 3)
   G4 = table(jc_angles.angles(4).KNEE_f,jc_angles.angles(4).KNEE_r,jc_angles.angles(4).KNEE_a,jc_angles.angles(4).ANKLE_f,jc_angles.angles(4).ANKLE_r,jc_angles.angles(4).ANKLE_a,jc_angles.angles(4).HIP_f,jc_angles.angles(4).HIP_r,jc_angles.angles(4).HIP_a);
end

writetable(G1,filename,'Sheet',1,'Range','B5','WriteVariableNames',false);

if(gait_cycle_count > 1)
    writetable(G2,filename,'Sheet',1,'Range','K5','WriteVariableNames',false);
end

if(gait_cycle_count > 2)
    writetable(G3,filename,'Sheet',1,'Range','T5','WriteVariableNames',false);
end

if(gait_cycle_count > 3)
    writetable(G4,filename,'Sheet',1,'Range','T5','WriteVariableNames',false);
end



end








