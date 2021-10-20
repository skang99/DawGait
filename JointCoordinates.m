function [BA_f,BA_r,BA_a,AP_f,AP_r,AP_a,BS_f,BS_r,BS_a,new_time,R_5th_M_z,frames,no_good_cycles,gait_cycle_count] = JointCoordinates(dynamic_trial,static_trial)

load(['Produced Data/' dynamic_trial])

gait_cycle_count = length(trial.gait_cycles);

no_good_cycles = 0;

R5M = trial.pos_data.R5M;
[R5M_x R5M_y R5M_z] = extract_XYZ(R5M);
R2M = trial.pos_data.R2M;
RLE = trial.pos_data.RLE;
RTR = trial.pos_data.RTR;
RME = trial.pos_data.RME;
RMS = trial.pos_data.RMS;
RGT = trial.pos_data.RGT;
RLO = trial.pos_data.RLO;
RLS = trial.pos_data.RLS;
T1 = trial.pos_data.T1;
Centroid = trial.pos_data.Centroid;
RDS = trial.pos_data.RDS;
RCR = trial.pos_data.RCR;
RAC = trial.pos_data.RAC;
ACB = trial.pos_data.ACB;
RSC1 = trial.pos_data.RSC1;
RCDS = trial.pos_data.RCDS;

%allows for static and dynamic matrix dimensions to agree
end_frame = length(trial.pos_data.R5M);

new_time = 1:end_frame;

[BA_f, BA_r, BA_a, BS_f, BS_r, BS_a, AP_f, AP_r, AP_a] = create_angle_data(RLE,RME,RGT,RLS,RMS,RLO,R5M,R2M,ACB,RCDS,RDS,RAC,RSC1);

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
    
    start_frame = trial.gait_cycles(i).start_frame;
    end_frame = trial.gait_cycles(i).end_frame - 1;
    gc_length = abs(start_frame - end_frame) + 1;
    
    if(trial.gait_cycles(i).RGT_RLE && trial.gait_cycles(i).RLO_RLS)
        %Elbow RGT_RLE & RLO_RLS
        disp("Creating elbow data for gc # " + i)
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
    
    %Carpus RLO_RLS & ACB_MP5
    if(trial.gait_cycles(i).ACB_R5M && trial.gait_cycles(i).RLO_RLS)
        disp("Creating carpus data for gc # " + i)
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
    
    %shoulder
    if(trial.gait_cycles(i).RGT_RLE && trial.gait_cycles(i).RSC1_RAC)
        disp("Creating shoulder data for gc # " + i)
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
       
    %all data included, excluding bad(segment check reject)/missing data
    x(i) = struct("ELB_f",tBA_f,"ELB_r",tBA_r,"ELB_a",tBA_a,"CARP_f",tAP_f,"CARP_r",tAP_r,"CARP_a",tAP_a,"SHLD_f",tBS_f,"SHLD_r",tBS_r,"SHLD_a",tBS_a,"start_frame",start_frame,"end_frame",end_frame);

end
frames = [];

for i = 1:length(trial.gait_cycles)
    frames(i) = trial.gait_cycles(i).start_frame;
    frames(i+1) = trial.gait_cycles(i).end_frame;
end

jc_angles = struct("trial_name",trial.trial_name,"angles",x,"ELB_f",nBA_f,"ELB_r",nBA_r,"ELB_a",nBA_a,"CARP_f",nAP_f,"CARP_r",nAP_r,"CARP_a",nAP_a,"SHLD_f",nBS_f,"SHLD_r",nBS_r,"SHLD_a",nBS_a,"normalized_angles",[]);

j = fieldnames(jc_angles.angles)';
j = j(1:9);

normalized_angles = zeros(9,100);

w = 1;
for i = 1:gait_cycle_count
    for k = 1:9
        tseries = jc_angles.angles(i).(j{k});
        if(tseries(1) ~= -1 && ~isnan(tseries(1)))
            yy = spline(1:length(tseries),tseries,linspace(1,length(tseries)));
            normalized_angles(w,:) = yy';
            w = w + 1;
        else
           normalized_angles(w,:) = (-1 * ones(100,1));
           w = w + 1;
        end
    end
end


temp_name = char(trial.trial_name + " Angles");

save(['Produced Data/' temp_name '.mat'],'jc_angles')

R_5th_M_z = R5M_z;

name = {trial.trial_name
        ''
        'GC#'};

filename = "Produced Data/" + trial.trial_name + " Angular Data.xlsx";

header = {'1' '' '' '' '' '' '' '' ''};
label = {'ELB_f' 'ELB_r' 'ELB_a' 'CARP_f' 'CARP_r' 'CARP_a' 'SHLD_f' 'SHLD_r' 'SHLD_a'};

header = repmat(header,1,gait_cycle_count);
label = repmat(label,1,gait_cycle_count);
index = 10;

for i=1:gait_cycle_count - 1
    header(index) = num2cell(i + 1);
    index = index + 9;
end

writecell(name,filename,'Range','A1');
writecell(header,filename,'Range','B3');
writecell(label,filename,'Range','B4');

G1 = table(normalized_angles');

writetable(G1,filename,'Sheet',1,'Range','B5','WriteVariableNames',false);

% G1 = table(jc_angles.angles(1).ELB_f,jc_angles.angles(1).ELB_r,jc_angles.angles(1).ELB_a,jc_angles.angles(1).CARP_f,jc_angles.angles(1).CARP_r,jc_angles.angles(1).CARP_a,jc_angles.angles(1).SHLD_f,jc_angles.angles(1).SHLD_r,jc_angles.angles(1).SHLD_a);
% 
% if(gait_cycle_count > 1)
%     G2 = table(jc_angles.angles(2).ELB_f,jc_angles.angles(2).ELB_r,jc_angles.angles(2).ELB_a,jc_angles.angles(2).CARP_f,jc_angles.angles(2).CARP_r,jc_angles.angles(2).CARP_a,jc_angles.angles(2).SHLD_f,jc_angles.angles(2).SHLD_r,jc_angles.angles(2).SHLD_a);
% end
% 
% if(gait_cycle_count > 2)
%     G3 = table(jc_angles.angles(3).ELB_f,jc_angles.angles(3).ELB_r,jc_angles.angles(3).ELB_a,jc_angles.angles(3).CARP_f,jc_angles.angles(3).CARP_r,jc_angles.angles(3).CARP_a,jc_angles.angles(3).SHLD_f,jc_angles.angles(3).SHLD_r,jc_angles.angles(3).SHLD_a);
% end
% 
% if(gait_cycle_count > 3)
%     G4 = table(jc_angles.angles(4).ELB_f,jc_angles.angles(4).ELB_r,jc_angles.angles(4).ELB_a,jc_angles.angles(4).CARP_f,jc_angles.angles(4).CARP_r,jc_angles.angles(4).CARP_a,jc_angles.angles(4).SHLD_f,jc_angles.angles(4).SHLD_r,jc_angles.angles(4).SHLD_a);
% end
% 
% writetable(G1,filename,'Sheet',1,'Range','B5','WriteVariableNames',false);
% 
% if(gait_cycle_count > 1)
%     writetable(G2,filename,'Sheet',1,'Range','K5','WriteVariableNames',false);
% end
% 
% if(gait_cycle_count > 2)
%     writetable(G3,filename,'Sheet',1,'Range','T5','WriteVariableNames',false);
% end
% 
% if(gait_cycle_count > 3)
%     writetable(G4,filename,'Sheet',1,'Range','AC5','WriteVariableNames',false);
% end


end








