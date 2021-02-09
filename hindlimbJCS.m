function [BA_f,BA_r,BA_a,AP_f,AP_r,AP_a,BS_f,BS_r,BS_a,new_time,R_5th_M_z,frames,no_good_cycles,gait_cycle_count] = hindlimbJCS(dynamic_trial,static_trial)

try
[sRMP5,sRGT,sRLEP,sRFH,sRLMA,sCRS,sRIWG,time,sRMEP,sRMMA,sRQUA,sGAS,sRCAL,sRISC] = hindlimbJCSextract(static_trial);
catch exception
    rethrow(exception)
end

load(['Produced Data/' dynamic_trial])

[static_5th_M_x static_5th_M_y static_5th_M_z] = extract_XYZ(sRMP5);
[static_RLEP_x static_RLEP_y static_RLEP_z] = extract_XYZ(sRLEP);
[static_Qua_x static_Qua_y static_Qua_z] = extract_XYZ(sRQUA);
[static_ME_x static_ME_y static_ME_z] = extract_XYZ(sRMEP);
[static_RGAS_x static_RGAS_y static_RGAS_z] = extract_XYZ(sGAS);
[static_RMMA_x static_RMMA_y static_RMMA_z] = extract_XYZ(sRMMA);
[static_RCAL_x static_RCAL_y static_RCAL_z] = extract_XYZ(sRCAL);
[static_CRS_x static_CRS_y static_CRS_z] = extract_XYZ(sCRS);
[static_RIWG_x static_RIWG_y static_RIWG_z] = extract_XYZ(sRIWG);
[static_RISC_x static_RISC_y static_RISC_z] = extract_XYZ(sRISC);

gait_cycle_count = length(trial.gait_cycles);

no_good_cycles = 0;

RMP5 = trial.pos_data.RMP5;
[R5M_x R5M_y R5M_z] = extract_XYZ(RMP5);

RMP2 = trial.pos_data.RMP2;
RLEP = trial.pos_data.RLEP;
RQUA = trial.pos_data.RQUA;
RMEP = trial.pos_data.RMEP;
RMMA = trial.pos_data.RMMA;
RGT = trial.pos_data.RGT;
RFH = trial.pos_data.RFH;
RLMA = trial.pos_data.RLMA;
CRS = trial.pos_data.CRS;
RIWG = trial.pos_data.RIWG;
RGAS = trial.pos_data.RGAS;
RCAL = trial.pos_data.RCAL;
RISC = trial.pos_data.RISC;

%allows for static and dynamic matrix dimensions to agree
end_frame = length(trial.pos_data.RMP5);

static_RLEP_x = mean(static_RLEP_x);
static_RLEP_y = mean(static_RLEP_y);
static_RLEP_z = mean(static_RLEP_y);

static_Qua_x = mean(static_Qua_x);
static_Qua_y = mean(static_Qua_y);
static_Qua_z = mean(static_Qua_z);

static_ME_x = mean(static_ME_x);
static_ME_y = mean(static_ME_y);
static_ME_z = mean(static_ME_z);

static_RGAS_x = mean(static_RGAS_x);
static_RGAS_y = mean(static_RGAS_y);
static_RGAS_z = mean(static_RGAS_z);

static_RMMA_x = mean(static_RMMA_x);
static_RMMA_y = mean(static_RMMA_y);
static_RMMA_z = mean(static_RMMA_z);

static_5th_M_x = mean(static_5th_M_x);
static_5th_M_y = mean(static_5th_M_y);
static_5th_M_z = mean(static_5th_M_z);

static_RCAL_x = mean(static_RCAL_x);
static_RCAL_y = mean(static_RCAL_y);
static_RCAL_z = mean(static_RCAL_z);

static_RISC_x = mean(static_RISC_x);
static_RISC_y = mean(static_RISC_y);
static_RISC_z = mean(static_RISC_z);

[RQUA_x RQUA_y RQUA_z] = extract_XYZ(RQUA);
[RGAS_x RGAS_y RGAS_z] = extract_XYZ(RGAS);
[RMMA_x RMMA_y RMMA_z] = extract_XYZ(RMMA);
[RMP2_x RMP2_y RMP2_z] = extract_XYZ(RMP2);

if(count_missing_frames(RMEP,20))
    disp("RMEP missing data")
    for k = 1:size(R5M_x)
        RMEP_x(k) = RQUA_x(k) - static_Qua_x + static_ME_x;
        RMEP_y(k) = RQUA_y(k) - static_Qua_y + static_ME_y;
        RMEP_z(k) = RQUA_z(k) - static_Qua_z + static_ME_z;
    end
end

if(count_missing_frames(RMMA,20))
    disp("RMMA missing data")
    for k = 1:size(R5M_x)
        RMMA_x(k) = RGAS_x(k) - static_RGAS_x + static_RMMA_x;
        RMMA_y(k) = RGAS_y(k) - static_RGAS_y + static_RMMA_y;
        RMMA_z(k) = RGAS_z(k) - static_RGAS_z + static_RMMA_z;
    end
end

if(count_missing_frames(RMP2,20))
    disp("RMP2 missing data")
    for k = 1:size(R5M_x)
        RMP2_x(k) = R5M_x(k) - static_5th_M_x + static_RCAL_x;
        RMP2_y(k) = R5M_y(k) - static_5th_M_y + static_RCAL_y;
        RMP2_z(k) = R5M_z(k) - static_5th_M_z + static_RCAL_z;  
    end
end

RMEP = [RMEP_x RMEP_y RMEP_z];
RMMA = [RMMA_x RMMA_y RMMA_z];
R2M = [RMP2_x RMP2_y RMP2_z];

new_time = 1:end_frame;

bm = [-1 -1 -1];

ps1 = 0;
ps2 = 0;
ps3 = 0;

[BA_f BA_r BA_a ps1 ps2 ps3] = create_angle_data(RLEP,RMEP,RGT,RLMA,RMMA,RFH,bm,ps1,ps2,ps3);
[AP_f AP_r AP_a ps1 ps2 ps3] = create_angle_data(bm,bm,bm,RMP5,R2M,RCAL,bm,ps1,ps2,ps3);
[BS_f BS_r BS_a ps1 ps2 ps3] = create_angle_data(RISC,CRS,RIWG,bm,bm,bm,bm,ps1,ps2,ps3);


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
    
    if(trial.gait_cycles(i).RGT_RLEP && trial.gait_cycles(i).RFH_RLMA)
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
    
    if(trial.gait_cycles(i).RCAL_RMP5 && trial.gait_cycles(i).RFH_RLMA)
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
    
    if(trial.gait_cycles(i).RGT_RLEP && trial.gait_cycles(i).RIWG_RISC && trial.gait_cycles(i).CRS_RISC)
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

R_5th_M_z = R5M_z;

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








