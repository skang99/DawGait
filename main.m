function [landmark_coord,trimmed_coord,cycle_time,gait_cycles] = main(dynamic_trial,name,static_trial,error_length,override_arr,show_graphs,gc_override_arr)

%Frontlimb trial function
[R5M,R2M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,time,RME,RMS,RTR,RCR,ACB,RAC,DLMC5] = create_frontlimb_data(dynamic_trial);

%Plots the position of the right 5th metacarpal during the dynamic trial
% figure(1) 
% plot(time,R_5th_M_z)
% xlabel('time(sec)')
% ylabel('position(m)')
% title('Z Coordinates of Paw in Time') 

cycle_time = time * 200;

R_5th_M_z = R5M(:,3);

landmark_coord = R_5th_M_z;

[cyc_start,cyc_end] = find_start_cycle_frame(R_5th_M_z);

R_5th_M_z = R_5th_M_z(cyc_start:cyc_end);


trimmed_coord = R_5th_M_z;

new_new_time = time(1:length(R_5th_M_z));

R5M = R5M(cyc_start:cyc_end,:);
R2M = R2M(cyc_start:cyc_end,:);
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
DLMC5 = DLMC5(cyc_start:cyc_end,:);

[gait_cycle_frame_locations, gait_cycle_count] = split_gait_cycles(R_5th_M_z);


if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
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

for n = 1:gait_cycle_count
   
    start_frame = gait_cycle_frame_locations(n);
    end_frame = gait_cycle_frame_locations(n+1);
    segment_error_checks = struct("RGT_RLE","RLO_RLS","RDS_Centroid","T1_Centroid","ACB_R5M",0);
    
    segment_error_checks.RGT_RLE = error_check(RGT,RLE,error_length,static_RGT_RLE,start_frame,end_frame,1);
    
    if(override_arr(oc))
        segment_error_checks.RGT_RLE = ~segment_error_checks.RGT_RLE;
    end
    
    segment_error_checks.RLO_RLS = error_check(RLO,RLS,error_length,static_RLO_RLS,start_frame,end_frame,0);
    
    if(override_arr(oc + 1))
        segment_error_checks.RLO_RLS = ~segment_error_checks.RLO_RLS;
    end
    

    segment_error_checks.RDS_Centroid = error_check(RDS,Centroid,error_length,static_RDS_Centroid,start_frame,end_frame,0);
    
    if(override_arr(oc + 2))
        segment_error_checks.RDS_Centroid = ~segment_error_checks.RDS_Centroid;
    end
    
    
    segment_error_checks.T1_Centroid = error_check(T1,Centroid,error_length,static_T1_Centroid,start_frame,end_frame,0);
    
    if(override_arr(oc + 3))    
        segment_error_checks.T1_Centroid = ~segment_error_checks.T1_Centroid;
    end
    
    
    segment_error_checks.ACB_R5M = error_check(ACB,R5M,error_length,static_ACB_R5M,start_frame,end_frame,0);
    
    if(override_arr(oc + 4))  
        segment_error_checks.ACB_R5M = ~segment_error_checks.ACB_R5M;
    end
    
    gait_cycles(n) = struct("start_frame",start_frame,"end_frame",end_frame,"seg_checks",segment_error_checks);
    
    oc = oc + 5;
end

if(show_graphs) 
    rgt_rle = lengthPlotter(static_RGT_RLE,RGT,RLE);

    new_new_time = 1:length(rgt_rle);

    %Plots the position of the right 5th metacarpal in one or more complete
    %gait cycles
    figure(1)
    plot(new_new_time,R_5th_M_z)
    xlabel('time(s)')
    ylabel('position(mm)')
    title('Gait Cycles')

    figure(2)
    plot(new_new_time,rgt_rle)
    hold on;
    yline(static_RGT_RLE);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RGT/RLE Static and Dynamic Segments vs Time')

    rlo_rls = lengthPlotter(static_RLO_RLS,RLO,RLS);
    figure(3)
    plot(new_new_time,rlo_rls)
    hold on;
    yline(static_RLO_RLS);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RLO/RLS Static and Dynamic Segments vs Time')

    rds_cent = lengthPlotter(static_RDS_Centroid,RDS,Centroid);
    figure(4)
    plot(new_new_time,rds_cent)
    hold on;
    yline(static_RDS_Centroid);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RDS/Cent Static and Dynamic Segments vs Time')

    t1_cent = lengthPlotter(static_T1_Centroid,T1,Centroid);
    figure(5)
    plot(new_new_time,t1_cent)
    hold on;
    yline(static_T1_Centroid);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of T1/Cent Static and Dynamic Segments vs Time')

    acb_r5m = lengthPlotter(static_ACB_R5M,ACB,R5M);
    figure(6)
    plot(new_new_time,acb_r5m)
    hold on;
    yline(static_ACB_R5M);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of ACB/R5M Static and Dynamic Segments vs Time')
end

for i = 1:gait_cycle_count
    disp("GC #" + i)
    disp(gait_cycles(i).seg_checks)
end

trial = struct("trial_name",-1,"gait_cycles",-1);

name = extractBetween(name,1,length(name)-4)
name = char(name);

trial.trial_name = name;

gc_override_array = gc_override_arr(1:gait_cycle_count);

for i = 1:length(gait_cycles)
    
    start_frame = gait_cycles(i).start_frame;
    end_frame = gait_cycles(i).end_frame;
    
    if(gc_override_array(i))
        gc_is_good = ~gait_cycles(i).seg_checks.RGT_RLE;
    else
      gc_is_good = gait_cycles(i).seg_checks.RGT_RLE;
    end
    
    tR5M = R5M(start_frame:end_frame,:);
    tR2M = R2M(start_frame:end_frame,:);
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
    tDLMC5 = DLMC5(start_frame:end_frame,:);
    
    x(i) = struct("is_good",gc_is_good,"R5M", tR5M, "R2M", tR2M, "RGT", tRGT, "RLE", tRLE, ...
    "RLO", tRLO, "RLS", tRLS, "T1", tT1, "RDS", tRDS, "Centroid", tCentroid, "RME", tRME, ...
    "RTR", tRTR, "RMS", tRMS, "DLMC5", tDLMC5, "start_frame", start_frame, "end_frame",end_frame, ...
    "RGT_RLE", gait_cycles(i).seg_checks.RGT_RLE, "RLO_RLS", gait_cycles(i).seg_checks.RLO_RLS, ...
    "RDS_Centroid", gait_cycles(i).seg_checks.RDS_Centroid, "T1_Centroid", gait_cycles(i).seg_checks.T1_Centroid, ...
    "ACB_R5M",gait_cycles(i).seg_checks.ACB_R5M);

end

trial.pos_data.R5M = R5M;
trial.pos_data.R2M = R2M;
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
trial.pos_data.DLMC5 = DLMC5;

trial.gait_cycles = x;

save(['Produced Data/' name '.mat'],'trial')

end
      








