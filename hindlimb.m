function [landmark_coord,trimmed_coord,cycle_time,gait_cycles] = hindlimb(dynamic_trial,name,static_trial,error_length,override_arr,show_graphs,gc_override_array)


[RMP5,RMP2,RGT,RLEP,RFH,RLMA,CRS,RIWG,time,RMEP,RMMA,RQUA,RGAS,RCAL,RISC] = create_hindlimb_data(dynamic_trial);

cycle_time = time * 200;

R_5th_M_z = RMP5(:,3);

landmark_coord = R_5th_M_z;

[cyc_start,cyc_end] = find_start_cycle_frame(R_5th_M_z);

R_5th_M_z = R_5th_M_z(cyc_start:cyc_end);
trimmed_coord = R_5th_M_z;

new_new_time = time(1:length(R_5th_M_z));

RMP5 = RMP5(cyc_start:cyc_end,:);
RMP2 = RMP2(cyc_start:cyc_end,:);
RGT = RGT(cyc_start:cyc_end,:);
RLEP = RLEP(cyc_start:cyc_end,:);
RFH = RFH(cyc_start:cyc_end,:);
RLMA = RLMA(cyc_start:cyc_end,:);
CRS = CRS(cyc_start:cyc_end,:);
RIWG = RIWG(cyc_start:cyc_end,:);
RMEP = RMEP(cyc_start:cyc_end,:);
RQUA = RQUA(cyc_start:cyc_end,:);
RMMA = RMMA(cyc_start:cyc_end,:); 
RGAS = RGAS(cyc_start:cyc_end,:);
RCAL = RCAL(cyc_start:cyc_end,:);
RISC = RISC(cyc_start:cyc_end,:);

[gait_cycle_frame_locations, gait_cycle_count] = split_gait_cycles(R_5th_M_z)

if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
end

if(gait_cycle_frame_locations == -1) 
    gait_cycle_frame_locations(1) = 1;
    gait_cycle_frame_locations(2) = cyc_end - cyc_start;
else
    gait_cycle_locations = gait_cycle_frame_locations;
end

try
    [static_RGT_RLEP,static_RFH_RLMA,static_RIWG_RISC,static_CRS_RISC,static_RCAL_RMP5] = create_static_hindlimb_data(static_trial,new_new_time,time);
catch exception  
    disp(getReport(exception))
    return
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
    segment_error_checks = struct("RGT_RLEP","RFH_RLMA","RIWG_RISC","CRS_RISC","RCAL_RMP5",0);
    
    segment_error_checks.RGT_RLEP = error_check(RGT,RLEP,error_length,static_RGT_RLEP,start_frame,end_frame,1);
    
    % RGT_RLEP
    if(override_arr(oc))
        segment_error_checks.RGT_RLEP = ~segment_error_checks.RGT_RLEP;
    end
    
    segment_error_checks.RFH_RLMA = error_check(RFH,RLMA,error_length,static_RFH_RLMA,start_frame,end_frame,0);
    
    % RFH_RLMA
    if(override_arr(oc + 1))
        segment_error_checks.RFH_RLMA = ~segment_error_checks.RFH_RLMA;
    end
    
    segment_error_checks.RIWG_RISC = error_check(RIWG,RISC,error_length,static_RIWG_RISC,start_frame,end_frame,0);
    
    % RIWG_RISC
    if(override_arr(oc + 2))
        segment_error_checks.RIWG_RISC = ~segment_error_checks.RIWG_RISC;
    end

    segment_error_checks.CRS_RISC = error_check(CRS,RISC,error_length,static_CRS_RISC,start_frame,end_frame,0);
    
    % CRS_RISC
    if(override_arr(oc + 3))    
        segment_error_checks.CRS_RISC = ~segment_error_checks.CRS_RISC;
    end
    
    segment_error_checks.RCAL_RMP5 = error_check(RCAL,RMP5,error_length,static_RCAL_RMP5,start_frame,end_frame,0);

    % RCAL_RMP5
    if(override_arr(oc + 4))  
        segment_error_checks.RCAL_RMP5 = ~segment_error_checks.RCAL_RMP5;
    end
    
    gait_cycles(n) = struct("start_frame",start_frame,"end_frame",end_frame,"seg_checks",segment_error_checks);
    
    oc = oc + 5;
 end

if(show_graphs) 
    rgt_rle = lengthPlotter(static_RGT_RLEP,RGT,RLEP);

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
    yline(static_RGT_RLEP);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RGT/RLEP Static and Dynamic Segments vs Time')

    rlo_rls = lengthPlotter(static_RFH_RLMA,RFH,RLMA);
    figure(3)
    plot(new_new_time,rlo_rls)
    hold on;
    yline(static_RFH_RLMA);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RFH/RLMA Static and Dynamic Segments vs Time')

    rds_cent = lengthPlotter(static_RIWG_RISC,RIWG,RISC);
    figure(4)
    plot(new_new_time,rds_cent)
    hold on;
    yline(static_RIWG_RISC);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RIWG/Cent Static and Dynamic Segments vs Time')

    t1_cent = lengthPlotter(static_CRS_RISC,CRS,RISC);
    figure(5)
    plot(new_new_time,t1_cent)
    hold on;
    yline(static_CRS_RISC);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of CRS/Cent Static and Dynamic Segments vs Time')

    acb_r5m = lengthPlotter(static_RCAL_RMP5,RCAL,RMP5);
    figure(6)
    plot(new_new_time,acb_r5m)
    hold on;
    yline(static_RCAL_RMP5);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of RCAL/RMP5 Static and Dynamic Segments vs Time')
end

for i = 1:gait_cycle_count
    disp("GC #" + i)
    disp(gait_cycles(i).seg_checks)
end

trial = struct("trial_name",-1,"gait_cycles",-1);

name = extractBetween(name,1,length(name)-4);
name = char(name);

trial.trial_name = name;


for i = 1:length(gait_cycles)
    
    start_frame = gait_cycles(i).start_frame;
    end_frame = gait_cycles(i).end_frame;
    
    if(gc_override_array(i))
        gc_is_good = ~gait_cycles(i).seg_checks.RGT_RLEP;
    else
      gc_is_good = gait_cycles(i).seg_checks.RGT_RLEP;
    end
    
    tRMP5 = RMP5(start_frame:end_frame,:);
    tRMP2 = RMP2(start_frame:end_frame,:);
    tRGT = RGT(start_frame:end_frame,:);
    tRLEP = RLEP(start_frame:end_frame,:);
    tRFH = RFH(start_frame:end_frame,:);
    tRLMA = RLMA(start_frame:end_frame,:);
    tCRS = CRS(start_frame:end_frame,:);
    tRIWG = RIWG(start_frame:end_frame,:);
    tRISC = RISC(start_frame:end_frame,:);
    tRMEP = RMEP(start_frame:end_frame,:); 
    tRQUA = RQUA(start_frame:end_frame,:);
    tRMMA = RMMA(start_frame:end_frame,:);
    
    x(i) = struct("is_good",gc_is_good,"RMP5", tRMP5, "RMP2", tRMP2,"RGT", tRGT, "RLEP", tRLEP, ...
    "RFH", tRFH, "RLMA", tRLMA, "CRS", tCRS, "RIWG", tRIWG, "RISC", tRISC, "RMEP", tRMEP, ...
    "RQUA", tRQUA, "RMMA", tRMMA,"start_frame", start_frame, "end_frame",end_frame, ...
    "RGT_RLEP", gait_cycles(i).seg_checks.RGT_RLEP, "RFH_RLMA", gait_cycles(i).seg_checks.RFH_RLMA, ...
    "RIWG_RISC", gait_cycles(i).seg_checks.RIWG_RISC, "CRS_RISC", gait_cycles(i).seg_checks.CRS_RISC, ...
    "RCAL_RMP5",gait_cycles(i).seg_checks.RCAL_RMP5);

    disp(x(i))

end

trial.pos_data.RMP5 = RMP5;
trial.pos_data.RMP2 = RMP2;
trial.pos_data.RGT = RGT;
trial.pos_data.RLEP = RLEP;
trial.pos_data.RFH = RFH;
trial.pos_data.RLMA = RLMA;
trial.pos_data.CRS = CRS;
trial.pos_data.RIWG = RIWG;
trial.pos_data.RISC = RISC;
trial.pos_data.RMEP = RMEP;
trial.pos_data.RQUA = RQUA;
trial.pos_data.RMMA = RMMA;
trial.pos_data.RGAS = RGAS;
trial.pos_data.RCAL = RCAL;
trial.pos_data.RISC = RISC;

trial.gait_cycles = x;

disp(name)

save(['Produced Data/' name '.mat'],'trial')

end
      








