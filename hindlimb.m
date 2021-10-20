function [landmark_coord,trimmed_coord,cycle_time,gait_cycles] = hindlimb(dynamic_trial,name,static_trial,error_length,override_arr,show_graphs,gc_override_array,reconstruct_arr,dir_adjust,show_reconst_graphs)

tic
[MP5,MP2,GT,LEP,FH,LMA,CRS,IWG,time,MEP,MMA,QUA,GAS,CAL,ISC,CGT,PTC,side] = create_hindlimb_data(dynamic_trial);
cycle_time = time * 200;

R_5th_M_z = MP5(:,3);
landmark_coord = R_5th_M_z;

[cyc_start,cyc_end] = find_start_cycle_frame(R_5th_M_z);

R_5th_M_z = R_5th_M_z(cyc_start:cyc_end);


MP5 = MP5(cyc_start:cyc_end,:); 
MP2 = MP2(cyc_start:cyc_end,:);
GT = GT(cyc_start:cyc_end,:);
LEP = LEP(cyc_start:cyc_end,:);
FH = FH(cyc_start:cyc_end,:);
LMA = LMA(cyc_start:cyc_end,:);
CRS = CRS(cyc_start:cyc_end,:);
IWG = IWG(cyc_start:cyc_end,:);
MEP = MEP(cyc_start:cyc_end,:);
QUA = QUA(cyc_start:cyc_end,:);
MMA = MMA(cyc_start:cyc_end,:); 
GAS = GAS(cyc_start:cyc_end,:);
CAL = CAL(cyc_start:cyc_end,:);
ISC = ISC(cyc_start:cyc_end,:);
CGT = CGT(cyc_start:cyc_end,:);
PTC = PTC(cyc_start:cyc_end,:);

[gait_cycle_frame_locations, gait_cycle_count] = split_gait_cycles(R_5th_M_z);

if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
end

[MP5_x, MP5_y, MP5_z] = extract_XYZ(MP5);

if(dir_adjust)
    if(any(diff(MP5_y) < 0))
        MP5(:,2) = MP5(:,2) * -1;
        MP2(:,2) = MP2(:,2) * -1;
        GT(:,2) = GT(:,2) * -1;
        LEP(:,2) = LEP(:,2) * -1;
        FH(:,2) = FH(:,2) * -1;
        LMA(:,2) = LMA(:,2) * -1;
        CRS(:,2) = CRS(:,2) * -1;
        IWG(:,2) = IWG(:,2) * -1;
        MEP(:,2) = MEP(:,2) * -1;
        QUA(:,2) = QUA(:,2) * -1;
        MMA(:,2) = MMA(:,2) * -1;
        GAS(:,2) = GAS(:,2) * -1;
        CAL(:,2) = CAL(:,2) * -1;
        ISC(:,2) = ISC(:,2) * -1;
        CGT(:,2) = CGT(:,2) * -1;
        PTC(:,2) = PTC(:,2) * -1;
    end
end

if(MP5_y(length(MP5_y)) > MP5_y(1))
    dir = 1;
else
    dir = -1;
end

if(gait_cycle_frame_locations == -1) 
    gait_cycle_frame_locations(1) = 1;
    gait_cycle_frame_locations(2) = cyc_end - cyc_start;
else
    gait_cycle_locations = gait_cycle_frame_locations;
end

new_new_time = time(1:length(R_5th_M_z));

try
    [static_GT_LEP,static_FH_LMA,static_IWG_ISC,static_CRS_ISC,static_CAL_MP5,sGT,sLEP,sMP5,sMP2,sFH,sLMA,sCRS,sIWG,sMEP,sQUA,sGAS,sISC,sMMA,sCAL,sDLMC5,sCGT,sPTC] = create_static_hindlimb_data(static_trial);
catch exception  
    disp(getReport(exception))
    return
end

gc = gait_cycle_locations;

MP5 = MP5(gc(1):gc(length(gc)),:);
MP2 = MP2(gc(1):gc(length(gc)),:);
GT = GT(gc(1):gc(length(gc)),:);
LEP = LEP(gc(1):gc(length(gc)),:);
FH = FH(gc(1):gc(length(gc)),:);
LMA = LMA(gc(1):gc(length(gc)),:);
CRS = CRS(gc(1):gc(length(gc)),:);
IWG = IWG(gc(1):gc(length(gc)),:);
MEP = MEP(gc(1):gc(length(gc)),:);
QUA = QUA(gc(1):gc(length(gc)),:);
MMA = MMA(gc(1):gc(length(gc)),:); 
GAS = GAS(gc(1):gc(length(gc)),:);
CAL = CAL(gc(1):gc(length(gc)),:);
ISC = ISC(gc(1):gc(length(gc)),:);
CGT = CGT(gc(1):gc(length(gc)),:);
PTC = PTC(gc(1):gc(length(gc)),:);

[QUA_x, QUA_y, QUA_z] = extract_XYZ(QUA);
[GT_x, GT_y, GT_z] = extract_XYZ(GT);
[LEP_x, LEP_y, LEP_z] = extract_XYZ(LEP);
[CGT_x, CGT_y, CGT_z] = extract_XYZ(CGT);
[MP2_x, MP2_y, MP2_z] = extract_XYZ(MP2);
[MP5_x, MP5_y, MP5_z] = extract_XYZ(MP5);
[GAS_x, GAS_y, GAS_z] = extract_XYZ(GAS);
[LMA_x, LMA_y, LMA_z] = extract_XYZ(LMA);
[FH_x, FH_y, FH_z] = extract_XYZ(FH);
[PTC_x, PTC_y, PTC_z] = extract_XYZ(PTC);
[CAL_x, CAL_y, CAL_z] = extract_XYZ(CAL);

new_new_time = time(1:length(MP5(:,3)));

%Readjusts GC locations to begin at frame #1 and end at cyc_end
if(gc(1) ~= 1)  
    gc = gc - gc(1);
    gc(1) = 1;
end

%Data for plotting GC graphs
trimmed_coord = MP5(:,3);

%MEP
if(reconstruct_arr(2))
    %Holds original MEP for later graphing
    org_MEP = MEP;
    MEP = [];
      
    [sQUA_x, sQUA_y, sQUA_z] = extract_XYZ(sQUA);
    [sGT_x, sGT_y, sGT_z] = extract_XYZ(sGT);
    [sLEP_x, sLEP_y, sLEP_z] = extract_XYZ(sLEP);
    [sMEP_x, sMEP_y, sMEP_z] = extract_XYZ(sMEP);
    [sCGT_x, sCGT_y, sCGT_z] = extract_XYZ(sCGT); 
    
    %Creates static segment lengths
    sQUA_x = mean(sQUA_x(1:100)); sQUA_y = mean(sQUA_y(1:100)); sQUA_z = mean(sQUA_z(1:100)); 
    sGT_x = mean(sGT_x(1:100)); sGT_y = mean(sGT_y(1:100)); sGT_z = mean(sGT_z(1:100));
    sLEP_x = mean(sLEP_x(1:100)); sLEP_y = mean(sLEP_y(1:100)); sLEP_z = mean(sLEP_z(1:100));
    sMEP_x = mean(sMEP_x(1:100)); sMEP_y = mean(sMEP_y(1:100)); sMEP_z = mean(sMEP_z(1:100));
    sCGT_x = mean(sCGT_x(1:100)); sCGT_y = mean(sCGT_y(1:100)); sCGT_z = mean(sCGT_z(1:100));
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(QUA(gc(i):gc(i+1)),GT(gc(i):gc(i+1)),LEP(gc(i):gc(i+1)),CGT(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            MEP = vertcat(MEP,zeros(gc(i+1)-gc(i),3));
            continue;
        end

        %If marker_good is 0, rm1_x,y,z will be 0
        %marker_count and marker_good should be consistent with one another
        marker_good = check_marker(QUA_x(gc(i):gc(i+1)));
        rm1_x = (QUA_x(gc(i):gc(i+1)) + -(dir * side)*abs(sQUA_x - sMEP_x)) * marker_good;
        rm1_y = (QUA_y(gc(i):gc(i+1)) + -(dir)*abs(sQUA_y - sMEP_y)) * marker_good;
        rm1_z = (QUA_z(gc(i):gc(i+1)) - abs(sQUA_z - sMEP_z)) * marker_good;
        
        marker_good = check_marker(GT_x(gc(i):gc(i+1)));
        rm2_x = (GT_x(gc(i):gc(i+1)) + -(dir * side)*abs(sGT_x - sMEP_x)) * marker_good;
        rm2_y = (GT_y(gc(i):gc(i+1)) + -(dir)*abs(sGT_y - sMEP_y)) * marker_good;
        rm2_z = (GT_z(gc(i):gc(i+1)) - abs(sGT_z - sMEP_z)) * marker_good;
        
        marker_good = check_marker(LEP_x(gc(i):gc(i+1)));
        rm3_x = (LEP_x(gc(i):gc(i+1)) + -(dir * side)*abs(sLEP_x - sMEP_x)) * marker_good;
        rm3_y = (LEP_y(gc(i):gc(i+1)) + -(dir)*abs(sLEP_y - sMEP_y)) * marker_good;
        rm3_z = (LEP_z(gc(i):gc(i+1)) - abs(sLEP_z - sMEP_z)) * marker_good;
        
        marker_good = check_marker(CGT_x(gc(i):gc(i+1)));
        rm4_x = (CGT_x(gc(i):gc(i+1)) + -(dir * side)*abs(sCGT_x - sMEP_x)) * marker_good;
        rm4_y = (CGT_y(gc(i):gc(i+1)) + -(dir)*abs(sCGT_y - sMEP_y)) * marker_good;
        rm4_z = (CGT_z(gc(i):gc(i+1)) - abs(sCGT_z - sMEP_z)) * marker_good;
        
        rMEP = [((rm1_x + rm2_x + rm3_x + rm4_x) / marker_count) ((rm1_y + rm2_y + rm3_y + rm4_y) / marker_count) ((rm1_z + rm2_z + rm3_z + rm4_z) / marker_count)];
        
        MEP = vertcat(MEP,rMEP);   
    end
    
    if(length(MP5) > length(MEP))
        MEP = vertcat(MEP,zeros(length(MP5)-length(MEP),3));
    else
        MEP = MEP(1:length(MP5),:);
    end
    
end

%MP2
if(reconstruct_arr(6))
    %Holds original MP2 for later graphing
    org_MP2 = MP2;
    MP2 = [];
  
    [sMP5_x, sMP5_y, sMP5_z] = extract_XYZ(sMP5);
    [sCAL_x, sCAL_y, sCAL_z] = extract_XYZ(sCAL);
    [sMP2_x, sMP2_y, sMP2_z] = extract_XYZ(sMP2);
 
    sMP5_x = mean(sMP5_x); sMP5_y = mean(sMP5_y); sMP5_z = mean(sMP5_z); 
    sCAL_x = mean(sCAL_x); sCAL_y = mean(sCAL_y); sCAL_z = mean(sCAL_z);
    sMP2_x = mean(sMP2_x); sMP2_y = mean(sMP2_y); sMP2_z = mean(sMP2_z);

    
  
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(MP5(gc(i):gc(i+1)),CAL(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            MP2 = vertcat(MP2,zeros(gc(i+1)-gc(i),3));
            continue;
        end

        marker_good = check_marker(MP5_x(gc(i):gc(i+1)));
        rm1_x = (MP5_x(gc(i):gc(i+1)) + -(dir * side)*abs(sMP5_x - sMP2_x)) * marker_good;
        rm1_y = (MP5_y(gc(i):gc(i+1)) + -(dir)*abs(sMP5_y - sMP2_y)) * marker_good;
        rm1_z = (MP5_z(gc(i):gc(i+1)) - abs(sMP5_z - sMP2_z)) * marker_good;
        
        marker_good = check_marker(CAL_x(gc(i):gc(i+1)));
        rm2_x = (CAL_x(gc(i):gc(i+1)) + -(dir * side)*abs(sCAL_x - sMP2_x)) * marker_good;
        rm2_y = (CAL_y(gc(i):gc(i+1)) + -(dir)*abs(sCAL_y - sMP2_y)) * marker_good;
        rm2_z = (CAL_z(gc(i):gc(i+1)) - abs(sCAL_z - sMP2_z)) * marker_good;
        
        
        rMP2 = [((rm1_x + rm2_x) / marker_count) ((rm1_y + rm2_y) / marker_count) ((rm1_z + rm2_z) / marker_count)];
        
        MP2 = vertcat(MP2,rMP2);   
    end
    
    if(length(MP5) > length(MP2))
        MP2 = vertcat(MP2,zeros(length(MP5)-length(MP2),3));
    else
        MP2 = MP2(1:length(MP5),:);
    end
end

%MMA
if(reconstruct_arr(4))  
    org_MMA = MMA;
    MMA = [];
     
    [sGAS_x, sGAS_y, sGAS_z] = extract_XYZ(sGAS);
    [sLMA_x, sLMA_y, sLMA_z] = extract_XYZ(sLMA);
    [sMMA_x, sMMA_y, sMMA_z] = extract_XYZ(sMMA);

    sGAS_x = mean(sGAS_x); sGAS_y = mean(sGAS_y); sGAS_z = mean(sGAS_z); 
    sLMA_x = mean(sLMA_x); sLMA_y = mean(sLMA_y); sLMA_z = mean(sLMA_z);
    sMMA_x = mean(sMMA_x); sMMA_y = mean(sMMA_y); sMMA_z = mean(sMMA_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(GAS(gc(i):gc(i+1)),LMA(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            MMA = vertcat(MMA,zeros(gc(i+1)-gc(i),3));
            continue;
        end

        marker_good = check_marker(GAS_x(gc(i):gc(i+1)));
        rm1_x = (GAS_x(gc(i):gc(i+1)) + -(dir * side)*abs(sGAS_x - sMMA_x)) * marker_good;
        rm1_y = (GAS_y(gc(i):gc(i+1)) + -(dir)*abs(sGAS_y - sMMA_y)) * marker_good;
        rm1_z = (GAS_z(gc(i):gc(i+1)) + abs(sGAS_z - sMMA_z)) * marker_good;

        marker_good = check_marker(LMA_x(gc(i):gc(i+1)));
        rm2_x = (LMA_x(gc(i):gc(i+1)) + -(dir * side)*abs(sLMA_x - sMMA_x)) * marker_good;
        rm2_y = (LMA_y(gc(i):gc(i+1)) + -(dir)*abs(sLMA_y - sMMA_y)) * marker_good;
        rm2_z = (LMA_z(gc(i):gc(i+1)) + abs(sLMA_z - sMMA_z)) * marker_good;

        rMMA = [((rm1_x + rm2_x) / marker_count) ((rm1_y + rm2_y) / marker_count) ((rm1_z + rm2_z) / marker_count)];
        
        MMA = vertcat(MMA,rMMA);   
    end
   
    
    if(length(MP5) > length(MMA))
        MMA = vertcat(MMA,zeros(length(MP5)-length(MMA),3));
    else
        MMA = MMA(1:length(MP5),:);
    end
    
    
end

%LEP
if(reconstruct_arr(1))   
    org_LEP = LEP;
    LEP = [];
   
    [sQUA_x, sQUA_y, sQUA_z] = extract_XYZ(sQUA);
    [sCGT_x, sCGT_y, sCGT_z] = extract_XYZ(sCGT);
    [sLEP_x, sLEP_y, sLEP_z] = extract_XYZ(sLEP);
    
    sQUA_x = mean(sQUA_x); sQUA_y = mean(sQUA_y); sQUA_z = mean(sQUA_z);
    sCGT_x = mean(sCGT_x); sCGT_y = mean(sCGT_y); sCGT_z = mean(sCGT_z);
    sLEP_x = mean(sLEP_x); sLEP_y = mean(sLEP_y); sLEP_z = mean(sLEP_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(QUA(gc(i):gc(i+1)),CGT(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            LEP = vertcat(LEP,zeros(gc(i+1)-gc(i),3));
            continue;
        end
               
        marker_good = check_marker(QUA_x(gc(i):gc(i+1)));
        rm1_x = (QUA_x(gc(i):gc(i+1)) + (dir * side)*abs(sQUA_x - sLEP_x)) * marker_good;
        rm1_y = (QUA_y(gc(i):gc(i+1)) + -(dir)*abs(sQUA_y - sLEP_y)) * marker_good;
        rm1_z = (QUA_z(gc(i):gc(i+1)) - abs(sQUA_z - sLEP_z)) * marker_good;
        
        marker_good = check_marker(CGT_x(gc(i):gc(i+1)));
        rm2_x = (CGT_x(gc(i):gc(i+1)) + (dir * side)*abs(sCGT_x - sLEP_x)) * marker_good;
        rm2_y = (CGT_y(gc(i):gc(i+1)) + -(dir)*abs(sCGT_y - sLEP_y)) * marker_good;
        rm2_z = (CGT_z(gc(i):gc(i+1)) - abs(sCGT_z - sLEP_z)) * marker_good;
       
        rLEP = [((rm1_x + rm2_x) / marker_count) ((rm1_y + rm2_y) / marker_count) ((rm1_z + rm2_z) / marker_count)];
        
        LEP = vertcat(LEP,rLEP);
    end
   
    if(length(MP5) > length(LEP))
        LEP = vertcat(LEP,zeros(length(MP5)-length(LEP),3));
    else
        LEP = LEP(1:length(MP5),:);
    end

end

%CAL
if(reconstruct_arr(3))    
    org_CAL = CAL;
    CAL = [];
    
    [sMP5_x, sMP5_y, sMP5_z] = extract_XYZ(sMP5);
    [sMP2_x, sMP2_y, sMP2_z] = extract_XYZ(sMP2);
    [sCAL_x, sCAL_y, sCAL_z] = extract_XYZ(sCAL);
       
    sMP5_x = mean(sMP5_x); sMP5_y = mean(sMP5_y); sMP5_z = mean(sMP5_z);
    sMP2_x = mean(sMP2_x); sMP2_y = mean(sMP2_y); sMP2_z = mean(sMP2_z);
    sCAL_x = mean(sCAL_x); sCAL_y = mean(sCAL_y); sCAL_z = mean(sCAL_z);
        
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(MP5(gc(i):gc(i+1)),MP2(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            CAL = vertcat(MP2,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        marker_good = check_marker(MP5_x(gc(i):gc(i+1)));
        rm1_x = (MP5_x(gc(i):gc(i+1)) + -(dir * side)*abs(sMP5_x - sCAL_x)) * marker_good;
        rm1_y = (MP5_y(gc(i):gc(i+1)) + -(dir)*abs(sMP5_y - sCAL_y)) * marker_good;
        rm1_z = (MP5_z(gc(i):gc(i+1)) - abs(sMP5_z - sCAL_z)) * marker_good;
        
        marker_good = check_marker(MP2_x(gc(i):gc(i+1)));
        rm2_x = (MP2_x(gc(i):gc(i+1)) + -(dir * side)*abs(sMP2_x - sCAL_x)) * marker_good;
        rm2_y = (MP2_y(gc(i):gc(i+1)) + -(dir)*abs(sMP2_y - sCAL_y)) * marker_good;
        rm2_z = (MP2_z(gc(i):gc(i+1)) - abs(sMP2_z - sCAL_z)) * marker_good;
        
        
        rMP2 = [((rm1_x + rm2_x) / marker_count) ((rm1_y + rm2_y) / marker_count) ((rm1_z + rm2_z) / marker_count)];
        
        CAL = vertcat(MP2,rMP2);
    end
    
    if(length(MP5) > length(CAL)) 
        CAL = vertcat(CAL,zeros(length(MP5)-length(CAL),3));
    else
        CAL = CAL(1:length(MP5),:);
    end
end

%FH
if(reconstruct_arr(5))
    org_FH = FH;
    FH = [];
    
    [PTC_x,PTC_y,PTC_z] = extract_XYZ(PTC);
    [sPTC_x, sPTC_y, sPTC_z] = extract_XYZ(sPTC);
    [sFH_x, sFH_y, sFH_z] = extract_XYZ(sFH);
   
    sPTC_x = mean(sPTC_x); sPTC_y = mean(sPTC_y); sPTC_z = mean(sPTC_z);
    sFH_x = mean(sFH_x); sFH_y = mean(sFH_y); sFH_z = mean(sFH_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(PTC(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            FH = vertcat(FH,zeros(gc(i+1)-gc(i),3));
            continue;
        end
                
        marker_good = check_marker(PTC_x(gc(i):gc(i+1)));
        rm1_x = (PTC_x(gc(i):gc(i+1)) + (dir * side)*abs(sPTC_x - sFH_x)) * marker_good;
        rm1_y = (PTC_y(gc(i):gc(i+1)) + -(dir)*abs(sPTC_y - sFH_y)) * marker_good;
        rm1_z = (PTC_z(gc(i):gc(i+1)) + abs(sPTC_z - sFH_z)) * marker_good;
           
        rFH = [rm1_x, rm1_y, rm1_z];
        
        FH = vertcat(FH,rFH);
    end
    
    if(length(MP5) > length(FH))
        FH = vertcat(FH,zeros(length(MP5)-length(FH),3));
    else
        FH = FH(1:length(MP5),:);
    end
      
end


if(show_reconst_graphs)
    %LEP
    if(reconstruct_arr(1))
        figure(1)
        plot(new_new_time,org_LEP(:,1))
        hold on;
        plot(new_new_time,QUA(:,1))
        plot(new_new_time,CGT(:,1))
        plot(new_new_time,LEP(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed LEP data')
        legend('original LEP','QUA','CGT','Reconstructed','Location','Northeast')
        
        figure(2)
        plot(new_new_time,org_LEP(:,2))
        hold on;
        plot(new_new_time,QUA(:,2))
        plot(new_new_time,CGT(:,2))
        plot(new_new_time,LEP(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed LEP data')
        legend('original LEP','QUA','CGT','Reconstructed','Location','Northeast')
        
        figure(3)
        plot(new_new_time,org_LEP(:,3))
        hold on;
        plot(new_new_time,QUA(:,3))
        plot(new_new_time,CGT(:,3))
        plot(new_new_time,LEP(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed LEP data')
        legend('original LEP','QUA','CGT','Reconstructed','Location','Northeast')
    end
    
    %MEP
    if(reconstruct_arr(2))
        figure(4)
        plot(new_new_time,org_MEP(:,1))
        hold on;
        plot(new_new_time,QUA(:,1))
        plot(new_new_time,GT(:,1))
        plot(new_new_time,LEP(:,1))
        plot(new_new_time,CGT(:,1))
        plot(new_new_time,MEP(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed MEP data')
        legend('original MEP','QUA','GT','LEP','CGT','Reconstructed','Location','Northeast')
        
        figure(5)
        plot(new_new_time,org_MEP(:,2))
        hold on;
        plot(new_new_time,QUA(:,2))
        plot(new_new_time,GT(:,2))
        plot(new_new_time,LEP(:,2))
        plot(new_new_time,CGT(:,2))
        plot(new_new_time,MEP(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed MEP data')
        legend('original MEP','QUA','GT','LEP','CGT','Reconstructed','Location','Northeast')
        
        figure(6)
        plot(new_new_time,org_MEP(:,3))
        hold on;
        plot(new_new_time,QUA(:,3))
        plot(new_new_time,GT(:,3))
        plot(new_new_time,LEP(:,3))
        plot(new_new_time,CGT(:,3))
        plot(new_new_time,MEP(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed MEP data')
        legend('original MEP','QUA','GT','LEP','CGT','Reconstructed','Location','Northeast')
    end
    
    %CAL
    if(reconstruct_arr(3))
        figure(7)
        plot(new_new_time,org_CAL(:,1))
        hold on;
        plot(new_new_time,MP5(:,1))
        plot(new_new_time,MP2(:,1))
        plot(new_new_time,CAL(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed CAL data')
        legend('original CAL','MP5','MP2','Reconstructed','Location','Northeast')
        
        figure(8)
        plot(new_new_time,org_CAL(:,2))
        hold on;
        plot(new_new_time,MP5(:,2))
        plot(new_new_time,MP2(:,2))
        plot(new_new_time,CAL(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed CAL data')
        legend('original CAL','MP5','MP2','Reconstructed','Location','Northeast')
        
        figure(9)
        plot(new_new_time,org_CAL(:,3))
        hold on;
        plot(new_new_time,MP5(:,3))
        plot(new_new_time,MP2(:,3))
        plot(new_new_time,CAL(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed CAL data')
        legend('original CAL','MP5','MP2','Reconstructed','Location','Northeast')
    end
    
    %FH
    if(reconstruct_arr(4))
        figure(10)
        plot(new_new_time,org_FH(:,1))
        hold on;
        plot(new_new_time,PTC(:,1))
        plot(new_new_time,FH(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed FH data')
        legend('original FH','PTC','Reconstructed','Location','Northeast')
        
        figure(11)
        plot(new_new_time,org_FH(:,2))
        hold on;
        plot(new_new_time,PTC(:,2))
        plot(new_new_time,FH(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed FH data')
        legend('original FH','PTC','Reconstructed','Location','Northeast')
        
        figure(12)
        plot(new_new_time,org_FH(:,3))
        hold on;
        plot(new_new_time,PTC(:,3))
        plot(new_new_time,FH(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed FH data')
        legend('original FH','PTC','Reconstructed','Location','Northeast')
    end
    
    %MMA
    if(reconstruct_arr(5))
        figure(13)
        plot(new_new_time,org_MMA(:,1))
        hold on;
        plot(new_new_time,GAS(:,1))
        plot(new_new_time,LMA(:,1))
        plot(new_new_time,MMA(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed ACB data')
        legend('original MMA','GAS','LMA','Reconstructed','Location','Northeast')
        
        figure(14)
        plot(new_new_time,org_MMA(:,2))
        hold on;
        plot(new_new_time,GAS(:,2))
        plot(new_new_time,LMA(:,2))
        plot(new_new_time,MMA(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed ACB data')
        legend('original MMA','GAS','LMA','Reconstructed','Location','Northeast')
        
        figure(15)
        plot(new_new_time,org_MMA(:,3))
        hold on;
        plot(new_new_time,GAS(:,3))
        plot(new_new_time,LMA(:,3))
        plot(new_new_time,MMA(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed ACB data')
        legend('original MMA','GAS','LMA','Reconstructed','Location','Northeast')
    end
end

  
% %Plots the position of the right 5th metacarpal during the static trial
% figure(3) %seperates graph into figure 1
% plot(time,static_R_5th_M_z) %plots original Z-paw data vs time
% xlabel('time(sec)') %label the x-axis
% ylabel('position(m)') %labels the y-axis
% title('Static Z Coordinates of Paw in Time') 

oc = 1;

for n = 1:gait_cycle_count
   
    start_frame = gc(n);
    end_frame = gc(n+1);
    
    segment_error_checks = struct("GT_LEP","FH_LMA","IWG_ISC",0,"CAL_MP5",0);
    
    segment_error_checks.GT_LEP = error_check(GT,LEP,error_length,static_GT_LEP,start_frame,end_frame,1);
    
    % GT_LEP
    if(override_arr(oc))
        segment_error_checks.GT_LEP = ~segment_error_checks.GT_LEP;
    end
    
    segment_error_checks.FH_LMA = error_check(FH,LMA,error_length,static_FH_LMA,start_frame,end_frame,0);
    
    % FH_LMA
    if(override_arr(oc + 1))
        segment_error_checks.FH_LMA = ~segment_error_checks.FH_LMA;
    end
    
    segment_error_checks.IWG_ISC = error_check(IWG,ISC,error_length,static_IWG_ISC,start_frame,end_frame,0);
    
    % IWG_ISC
    if(override_arr(oc + 2))
        segment_error_checks.IWG_ISC = ~segment_error_checks.IWG_ISC;
    end
   
    segment_error_checks.CAL_MP5 = error_check(CAL,MP5,error_length,static_CAL_MP5,start_frame,end_frame,0);

    % CAL_MP5
    if(override_arr(oc + 3))  
        segment_error_checks.CAL_MP5 = ~segment_error_checks.CAL_MP5;
    end
    
    gait_cycles(n) = struct("start_frame",start_frame,"end_frame",end_frame,"seg_checks",segment_error_checks);
    
    oc = oc + 4;
end

if(show_graphs) 
    GT_rle = lengthPlotter(static_GT_LEP,GT,LEP);

    new_new_time = 1:length(GT_rle);
    
    figure(17)
    plot(new_new_time,GT_rle)
    hold on;
    yline(static_GT_LEP);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Femur Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')

    fh_rlma = lengthPlotter(static_FH_LMA,FH,LMA);
    figure(18)
    plot(new_new_time,fh_rlma)
    hold on;
    yline(static_FH_LMA);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Tibular/Fibular Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')

    rds_cent = lengthPlotter(static_IWG_ISC,IWG,ISC);
    figure(19)
    plot(new_new_time,rds_cent)
    hold on;
    yline(static_IWG_ISC);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Pelvis Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')

    acb_r5m = lengthPlotter(static_CAL_MP5,CAL,MP5);
    figure(21)
    plot(new_new_time,acb_r5m)
    hold on;
    yline(static_CAL_MP5);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Tarsus Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')
end

% for i = 1:gait_cycle_count
%     disp("GC #" + i)
%     disp(gait_cycles(i).seg_checks)
% end

trial = struct("trial_name",-1,"gait_cycles",-1);

name = extractBetween(name,1,length(name)-4);
name = char(name);

if(side == -1)
    trial.side = "left";
else
    trial.side = "right";
end %if 

trial.trial_name = name;


for i = 1:length(gait_cycles)
    
    start_frame = gait_cycles(i).start_frame;
    end_frame = gait_cycles(i).end_frame;
    
    if(gc_override_array(i))
        gc_is_good = ~gait_cycles(i).seg_checks.GT_LEP;
    else
      gc_is_good = gait_cycles(i).seg_checks.GT_LEP;
    end
    
    tMP5 = MP5(start_frame:end_frame,:);
    tMP2 = MP2(start_frame:end_frame,:);
    tGT = GT(start_frame:end_frame,:);
    tLEP = LEP(start_frame:end_frame,:);
    tFH = FH(start_frame:end_frame,:);
    tLMA = LMA(start_frame:end_frame,:);
    tCRS = CRS(start_frame:end_frame,:);
    tIWG = IWG(start_frame:end_frame,:);
    tISC = ISC(start_frame:end_frame,:);
    tMEP = MEP(start_frame:end_frame,:); 
    tQUA = QUA(start_frame:end_frame,:);
    tMMA = MMA(start_frame:end_frame,:);
    
    x(i) = struct("is_good",gc_is_good,"MP5", tMP5, "MP2", tMP2,"GT", tGT, "LEP", tLEP, ...
    "FH", tFH, "LMA", tLMA, "CRS", tCRS, "IWG", tIWG, "ISC", tISC, "MEP", tMEP, ...
    "QUA", tQUA, "MMA", tMMA,"start_frame", start_frame, "end_frame",end_frame, ...
    "GT_LEP", gait_cycles(i).seg_checks.GT_LEP, "FH_LMA", gait_cycles(i).seg_checks.FH_LMA, ...
    "IWG_ISC", gait_cycles(i).seg_checks.IWG_ISC, ...
    "CAL_MP5",gait_cycles(i).seg_checks.CAL_MP5);
end

trial.pos_data.MP5 = MP5;
trial.pos_data.MP2 = MP2;
trial.pos_data.GT = GT;
trial.pos_data.LEP = LEP;
trial.pos_data.FH = FH;
trial.pos_data.LMA = LMA;
trial.pos_data.CRS = CRS;
trial.pos_data.IWG = IWG;
trial.pos_data.ISC = ISC;
trial.pos_data.MEP = MEP;
trial.pos_data.QUA = QUA;
trial.pos_data.MMA = MMA;
trial.pos_data.GAS = GAS;
trial.pos_data.CAL = CAL;
trial.pos_data.ISC = ISC;

trial.gait_cycles = x;

save(['Produced Data/' name '.mat'],'trial')
toc

end
      








