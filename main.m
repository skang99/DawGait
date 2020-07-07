function [landmark_coord,trimmed_coord,cycle_time,gait_cycles] = main(dynamic_trial,static_trial,error_length,override_arr)

% try
%     [R5M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,RME,RTR,RMS,time] = create_gait_cycles(dynamic_trial);
% catch exception %#ok<NASGU>;
%     return;
% end

[R5M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,time,RME,RMS,RTR,RCR,ACB,RAC] = create_gait_cycles(dynamic_trial);

%Plots the position of the right 5th metacarpal during the dynamic trial
% figure(1) 
% plot(time,R_5th_M_z)
% xlabel('time(sec)')
% ylabel('position(m)')
% title('Z Coordinates of Paw in Time') 

cycle_time = time;

R_5th_M_z = R5M(:,3);

landmark_coord = R_5th_M_z;

[cyc_start,cyc_end] = find_start_cycle_frame(R_5th_M_z);

R_5th_M_z = R_5th_M_z(cyc_start:cyc_end);
trimmed_coord = R_5th_M_z;

new_new_time = time(1:length(R_5th_M_z));

%Plots the position of the right 5th metacarpal in one or more complete
%gait cycles
figure(1)
plot(new_new_time,R_5th_M_z)
xlabel('time(s)')
ylabel('position(mm)')
title('Gait Cycles')

R5M = R5M(cyc_start:cyc_end,:);
RGT = RGT(cyc_start:cyc_end,:);
RLE = RLE(cyc_start:cyc_end,:);
RLO = RLO(cyc_start:cyc_end,:);
RLS = RLS(cyc_start:cyc_end,:);
T1 = T1(cyc_start:cyc_end,:);
RDS = RDS(cyc_start:cyc_end,:);
Centroid = Centroid(cyc_start:cyc_end,:);
RME = RME(cyc_start:cyc_end,:);
RTR = RTR(cyc_start:cyc_end,:);
RMS = RMS(cyc_start:cyc_end,:); 
RCR = RCR(cyc_start:cyc_end,:);
ACB = ACB(cyc_start:cyc_end,:);
RAC = RAC(cyc_start:cyc_end,:);

[gait_cycle_frame_locations, gait_cycle_count] = find_gait_cycle_frames(R_5th_M_z);


gait_cycle_count = gait_cycle_count + 1;

if(gait_cycle_frame_locations == -1) 
    gait_cycle_frame_locations(1) = 1;
    gait_cycle_frame_locations(2) = cyc_end - cyc_start;
else
    gait_cycle_locations = gait_cycle_frame_locations;
end

try
    [static_RGT_RLE,static_RLO_RLS,static_RDS_Centroid,static_T1_Centroid,static_ACB_R5M] = create_static_data(static_trial,new_new_time,time);
catch exception  
    disp(getReport(exception))
    
end
                   
% %Plots the position of the right 5th metacarpal during the static trial
% figure(3) %seperates graph into figure 1
% plot(time,static_R_5th_M_z) %plots original Z-paw data vs time
% xlabel('time(sec)') %label the x-axis
% ylabel('position(m)') %labels the y-axis
% title('Static Z Coordinates of Paw in Time') 

oc = 1;


for n = 1:gait_cycle_count - 1
   
    start_frame = gait_cycle_frame_locations(n);
    end_frame = gait_cycle_frame_locations(n+1);
    segment_error_checks = struct("RGT_RLE","RLO_RLS","RDS_Centroid","T1_Centroid","ACB_R5M",0);
    
    if(override_arr(oc))
        segment_error_checks.RGT_RLE = 1;
    else
        segment_error_checks.RGT_RLE = error_check(RGT,RLE,error_length,static_RGT_RLE,start_frame,end_frame,0);
    end
    

    if(override_arr(oc + 1))
        segment_error_checks.RLO_RLS = 1;
    else
        segment_error_checks.RLO_RLS = error_check(RLO,RLS,error_length,static_RLO_RLS,start_frame,end_frame,0);
    end
    

    if(override_arr(oc + 2))
        segment_error_checks.RDS_Centroid = 1;
    else
        segment_error_checks.RDS_Centroid = error_check(RDS,Centroid,error_length,static_RDS_Centroid,start_frame,end_frame,0);
    end
    

    if(override_arr(oc + 3))    
        segment_error_checks.T1_Centroid = 1;
    else
     segment_error_checks.T1_Centroid = error_check(T1,Centroid,error_length,static_T1_Centroid,start_frame,end_frame,0);
    end
    

    if(override_arr(oc + 4))  
        segment_error_checks.ACB_R5M = 1;
    else
        segment_error_checks.ACB_R5M = error_check(ACB,R5M,error_length,static_ACB_R5M,start_frame,end_frame,0);
    end
    
    gait_cycles(n) = struct("start_frame",start_frame,"end_frame",end_frame,"seg_checks",segment_error_checks);
    
    oc = oc + 5;
end


rgt_rle = lengthPlotter(static_RGT_RLE,RGT,RLE);

new_new_time = 1:length(rgt_rle);

figure(2)
plot(new_new_time,rgt_rle)
hold on;
yline(static_RGT_RLE);
hold off;
xlabel('frames')
ylabel('length(mm)')
title('Length of RGT/RLE Static and Dynamic Segments vs Time')

rlo_rls = lengthPlotter(static_RLO_RLS,RLO,RLS);
figure(2)
plot(new_new_time,rlo_rls)
hold on;
yline(static_RLO_RLS);
hold off;
xlabel('frames')
ylabel('length(mm)')
title('Length of RLO/RLS Static and Dynamic Segments vs Time')

rds_cent = lengthPlotter(static_RDS_Centroid,RDS,Centroid);
figure(3)
plot(new_new_time,rds_cent)
hold on;
yline(static_RDS_Centroid);
hold off;
xlabel('frames')
ylabel('length(mm)')
title('Length of RDS/Cent Static and Dynamic Segments vs Time')

t1_cent = lengthPlotter(static_T1_Centroid,T1,Centroid);
figure(4)
plot(new_new_time,t1_cent)
hold on;
yline(static_T1_Centroid);
hold off;
xlabel('frames')
ylabel('length(mm)')
title('Length of T1/Cent Static and Dynamic Segments vs Time')

acb_r5m = lengthPlotter(static_ACB_R5M,ACB,R5M);
figure(5)
plot(new_new_time,acb_r5m)
hold on;
yline(static_ACB_R5M);
hold off;
xlabel('frames')
ylabel('length(mm)')
title('Length of ACB/R5M Static and Dynamic Segments vs Time')

for i = 1:gait_cycle_count - 1
    disp("GC #" + i)
    disp(gait_cycles(i).seg_checks)
end

trial = struct("trial_name",-1,"gait_cycles",-1);

%temp fix, not all names are the same length, should parse for the first
%occurence of a number
name = extractBetween(dynamic_trial,8,strlength(dynamic_trial)-4);
trial.trial_name = name;

for i = 1:length(gait_cycles)
    
    start_frame = gait_cycles(i).start_frame;
    end_frame = gait_cycles(i).end_frame;
    gc_is_good = gait_cycles(i).seg_checks.RGT_RLE;
    
    tR5M = R5M(start_frame:end_frame,:);
    tRGT = RGT(start_frame:end_frame,:);
    tRLE = RLE(start_frame:end_frame,:);
    tRLO = RLO(start_frame:end_frame,:);
    tRLS = RLS(start_frame:end_frame,:);
    tT1 = T1(start_frame:end_frame,:);
    tRDS = RDS(start_frame:end_frame,:);
    tCentroid = Centroid(start_frame:end_frame,:);
    tRME = RME(start_frame:end_frame,:); 
    tRTR = RTR(start_frame:end_frame,:);
    tRMS = RMS(start_frame:end_frame,:);
    
    x(i) = struct("is_good",gc_is_good,"R5M", tR5M, "RGT", tRGT, "RLE", tRLE, ...
    "RLO", tRLO, "RLS", tRLS, "T1", tT1, "RDS", tRDS, "Centroid", tCentroid, "RME", tRME, ...
    "RTR", tRTR, "RMS", tRMS, "start_frame", start_frame, "end_frame",end_frame, ...
    "RGT_RLE", gait_cycles(i).seg_checks.RGT_RLE, "RLO_RLS", gait_cycles(i).seg_checks.RLO_RLS, ...
    "RDS_Centroid", gait_cycles(i).seg_checks.RDS_Centroid, "T1_Centroid", gait_cycles(i).seg_checks.T1_Centroid, ...
    "ACB_R5M",gait_cycles(i).seg_checks.ACB_R5M);

end

trial.pos_data.R5M = R5M;
trial.pos_data.RGT = RGT;
trial.pos_data.RLE = RLE;
trial.pos_data.RLO = RLO;
trial.pos_data.RLS = RLS;
trial.pos_data.T1 = T1;
trial.pos_data.RDS = RDS;
trial.pos_data.Centroid = Centroid;
trial.pos_data.RME = RME;
trial.pos_data.RTR = RTR;
trial.pos_data.RMS = RMS;
trial.pos_data.RCR = RCR;
trial.pos_data.ACB = ACB;
trial.pos_data.RAC = RAC;

trial.gait_cycles = x;

name = char(name);

save([name '.mat'],'trial')

end
      








