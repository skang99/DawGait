function [landmark_coord,trimmed_coord,cycle_time,gait_cycles] = main(dynamic_trial,name,static_trial,error_length,override_arr,show_graphs,gc_override_arr,reconstruct_arr,dir_adjust,show_reconst_graphs)

tic
%Frontlimb trial function
[R5M,R2M,RGT,RLE,RLO,RLS,T1,RDS,Centroid,time,RME,RMS,RTR,RCR,ACB,RAC,DLMC5,VTR1,VTR2,VTR3,RSC1,RCDS] = create_frontlimb_data(dynamic_trial);

cycle_time = time * 200;

R_5th_M_z = R5M(:,3);

landmark_coord = R_5th_M_z;

[cyc_start,cyc_end] = find_start_cycle_frame(R_5th_M_z);

R_5th_M_z = R_5th_M_z(cyc_start:cyc_end);

trimmed_coord = R_5th_M_z;

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
VTR1 = VTR1(cyc_start:cyc_end,:);
VTR2 = VTR2(cyc_start:cyc_end,:);
VTR3 = VTR3(cyc_start:cyc_end,:);
RCDS = RCDS(cyc_start:cyc_end,:);
RSC1 = RSC1(cyc_start:cyc_end,:);

[gait_cycle_frame_locations, gait_cycle_count] = split_gait_cycles(R_5th_M_z);
gc = gait_cycle_frame_locations;

if(gait_cycle_count == 0)
    disp("No gait cycles found")
    return
end

try
    [static_RGT_RLE,static_RLO_RLS,static_RSC1_RAC,static_RDS_RAC,static_ACB_R5M,sT1,sVTR1,sSC1,sDLMC5,sACB,sRDS,sRTR,sRME,sRCR,sR5M,sR2M,sRMS,sRGT,sRLE,sCentroid,sRLS,sRLO,sVTR2,sVTR3,sRCDS] = create_static_data(static_trial);
catch exception  
    disp(getReport(exception))
end

R5M = R5M(gc(1):gc(length(gc)),:);
R2M = R2M(gc(1):gc(length(gc)),:);
RGT = RGT(gc(1):gc(length(gc)),:);
RLE = RLE(gc(1):gc(length(gc)),:);
RLO = RLO(gc(1):gc(length(gc)),:);
RLS = RLS(gc(1):gc(length(gc)),:);
T1 = T1(gc(1):gc(length(gc)),:);
RDS = RDS(gc(1):gc(length(gc)),:);
Centroid = Centroid(gc(1):gc(length(gc)),:);
RME = RME(gc(1):gc(length(gc)),:);
RTR = RTR(gc(1):gc(length(gc)),:);
RMS = RMS(gc(1):gc(length(gc)),:); 
RCR = RCR(gc(1):gc(length(gc)),:);
ACB = ACB(gc(1):gc(length(gc)),:);
RAC = RAC(gc(1):gc(length(gc)),:);
DLMC5 = DLMC5(gc(1):gc(length(gc)),:);
VTR1 = VTR1(gc(1):gc(length(gc)),:);
VTR2 = VTR2(gc(1):gc(length(gc)),:);
VTR3 = VTR3(gc(1):gc(length(gc)),:);
RCDS = RCDS(gc(1):gc(length(gc)),:);
RSC1 = RSC1(gc(1):gc(length(gc)),:);

new_new_time = time(1:length(R5M(:,3)));

[RTR_x, RTR_y, RTR_z] = extract_XYZ(RTR);
[RLS_x, RLS_y, RLS_z] = extract_XYZ(RLS);
[RLO_x, RLO_y, RLO_z] = extract_XYZ(RLO);
[RGT_x, RGT_y, RGT_z] = extract_XYZ(RGT);
[RLE_x, RLE_y, RLE_z] = extract_XYZ(RLE);
[RMS_x, RMS_y, RMS_z] = extract_XYZ(RMS);
[RCR_x, RCR_y, RCR_z] = extract_XYZ(RCR);
[R5M_x, R5M_y, R5M_z] = extract_XYZ(R5M);
[R2M_x, R2M_y, R2M_z] = extract_XYZ(R2M);
[T1_x, T1_y, T1_z] = extract_XYZ(T1);
[RME_x, RME_y, RME_z] = extract_XYZ(RME);
[ACB_x, ACB_y, ACB_z] = extract_XYZ(ACB);
[VTR1_x, VTR1_y, VTR1_z] = extract_XYZ(VTR1);
[Centroid_x, Centroid_y, Centroid_z] = extract_XYZ(Centroid);
[VTR1_x, VTR1_y, VTR1_z] = extract_XYZ(VTR1);
[VTR2_x, VTR2_y, VTR2_z] = extract_XYZ(VTR2);

if(isempty(VTR3))
    VTR3 = zeros(length(R5M_x),3);   
    sVTR3 = zeros(length(R5M_x),3);
end

[VTR3_x, VTR3_y, VTR3_z] = extract_XYZ(VTR3);
[RCDS_x, RCDS_y, RCDS_z] = extract_XYZ(RCDS);

%Readjusts GC locations to begin at frame #1 and end at cyc_end
if(gc(1) ~= 1)  
    gc = gc - gc(1);
    gc(1) = 1;
end

%Data for plotting GC graphs
trimmed_coord = R5M_z;

if(isempty(DLMC5))
    DLMC5 = zeros(length(R5M_x),3);
    sDLMC5 = zeros(length(R5M_x),3);
end

[DLMC5_x, DLMC5_y, DLMC5_z] = extract_XYZ(DLMC5);

%Automatic direction adjustment
if(dir_adjust)
    if(any(diff(R5M_y) < 0))
        R5M(:,2) = R5M(:,2) * -1;
        R2M(:,2) = R2M(:,2) * -1;
        RGT(:,2) = RGT(:,2) * -1;
        RLE(:,2) = RLE(:,2) * -1;
        RLO(:,2) = RLO(:,2) * -1;
        RLS(:,2) = RLS(:,2) * -1;
        T1(:,2) = T1(:,2) * -1;
        RDS(:,2) = RDS(:,2) * -1;
        Centroid(:,2) = Centroid(:,2) * -1;
        RME(:,2) = RME(:,2) * -1;
        RTR(:,2) = RTR(:,2) * -1;
        RMS(:,2) = RMS(:,2) * -1;
        RCR(:,2) = RCR(:,2) * -1;
        ACB(:,2) = ACB(:,2) * -1;
        RAC(:,2) = RAC(:,2) * -1;
        DLMC5(:,2) = DLMC5(:,2) * -1;
    end
end

if(R5M_y(length(R5M_y)) > R5M_y(1))
    dir = 1;
else
    dir = -1;
end

%ACB
if(reconstruct_arr(5))  
    org_ACB = ACB;
    ACB = [];
    
    [sR2M_x, sR2M_y, sR2M_z] = extract_XYZ(sR2M);
    [sR5M_x, sR5M_y, sR5M_z] = extract_XYZ(sR5M);
    [sDLMC5_x, sDLMC5_y, sDLMC5_z] = extract_XYZ(sDLMC5);    
    [sACB_x, sACB_y, sACB_z] = extract_XYZ(sACB);
    
    sR2M_x = mean(sR2M_x); sR2M_y = mean(sR2M_y); sR2M_z = mean(sR2M_z); 
    sR5M_x = mean(sR5M_x); sR5M_y = mean(sR5M_y); sR5M_z = mean(sR5M_z); 
    sDLMC5_x = mean(sDLMC5_x); sDLMC5_y = mean(sDLMC5_y); sDLMC5_z = mean(sDLMC5_z);     
    sACB_x = mean(sACB_x); sACB_y = mean(sACB_y); sACB_z = mean(sACB_z); 
   
    for i = 1:gait_cycle_count        
        marker_count = count_good_markers(R2M_x(gc(i):gc(i+1)),R5M_x(gc(i):gc(i+1)),DLMC5_x(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            ACB = vertcat(ACB,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        marker_good = check_marker(R2M_x(gc(i):gc(i+1)));        
        rm1_x = (R2M_x(gc(i):gc(i+1))+ -(dir)*abs(sR2M_x - sACB_x)) * marker_good;
        rm1_y = (R2M_y(gc(i):gc(i+1)) + -(dir)*abs(sR2M_y - sACB_y)) * marker_good;
        rm1_z = (R2M_z(gc(i):gc(i+1)) - abs(sR2M_z - sACB_z)) * marker_good;
        
        marker_good = check_marker(R5M_x(gc(i):gc(i+1)));    
        rm2_x = (R5M_x(gc(i):gc(i+1))+ -(dir)*abs(sR5M_x - sACB_x)) * marker_good;
        rm2_y = (R5M_y(gc(i):gc(i+1)) + -(dir)*abs(sR5M_y - sACB_y)) * marker_good;
        rm2_z = (R5M_z(gc(i):gc(i+1)) - abs(sR5M_z - sACB_z)) * marker_good;
        
        marker_good = check_marker(DLMC5_x(gc(i):gc(i+1)));    
        rm3_x = (DLMC5_x(gc(i):gc(i+1))+ -(dir)*abs(sDLMC5_x - sACB_x)) * marker_good;
        rm3_y = (DLMC5_y(gc(i):gc(i+1)) + -(dir)*abs(sDLMC5_y - sACB_y)) * marker_good;
        rm3_z = (DLMC5_z(gc(i):gc(i+1)) - abs(sDLMC5_z - sACB_z)) * marker_good;
        
      
        rACB = [((rm1_x + rm2_x + rm3_x) / marker_count) ((rm1_y + rm2_y + rm3_y) / marker_count) ((rm1_z + rm2_z + rm3_z) / marker_count)];
        
        ACB = vertcat(ACB,rACB);
        
    end
   
    if(length(R5M) > length(ACB))
        ACB = vertcat(ACB,zeros(length(R5M)-length(ACB),3));
    else
        ACB = ACB(1:length(R5M),:);
    end
end

%RDS
if(reconstruct_arr(3))
    org_RDS = RDS;
    RDS = [];
    
    [sCentroid_x, sCentroid_y, sCentroid_z] = extract_XYZ(sCentroid);
    sCentroid_x = mean(sCentroid_x); sCentroid_y = mean(sCentroid_y); sCentroid_z = mean(sCentroid_z); 
    
    [sRCDS_x, sRCDS_y, sRCDS_z] = extract_XYZ(sRCDS);
    sRCDS_x = mean(sRCDS_x); sRCDS_y = mean(sRCDS_y); sRCDS_z = mean(sRCDS_z); 

    [sRDS_x, sRDS_y, sRDS_z] = extract_XYZ(sRDS);
    sRDS_x = mean(sRDS_x); sRDS_y = mean(sRDS_y); sRDS_z = mean(sRDS_z); 
    
    for i = 1:gait_cycle_count      
        marker_count = count_good_markers(Centroid(gc(i):gc(i+1)),RCDS(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            RDS = vertcat(RDS,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        marker_good = check_marker(Centroid_x(gc(i):gc(i+1)));
        rm1_x = (Centroid_x(gc(i):gc(i+1))+ -(dir)*abs(sCentroid_x - sRDS_x)) * marker_good;
        rm1_y = (Centroid_y(gc(i):gc(i+1)) + -(dir)*abs(sCentroid_y - sRDS_y)) * marker_good;
        rm1_z = (Centroid_z(gc(i):gc(i+1)) - abs(sCentroid_z - sRDS_z)) * marker_good;
        
        marker_good = check_marker(RCDS_x(gc(i):gc(i+1)));
        rm2_x = (RCDS_x(gc(i):gc(i+1)) + -(dir)*abs(sRCDS_x - sRDS_x)) * marker_good;
        rm2_y = (RCDS_y(gc(i):gc(i+1)) + -(dir)*abs(sRCDS_y - sRDS_y)) * marker_good;
        rm2_z = (RCDS_z(gc(i):gc(i+1)) + abs(sRCDS_z - sRDS_z)) * marker_good;
        
        rRDS = [((rm1_x + rm2_x) / marker_count) ((rm1_y + rm2_y) / marker_count) ((rm1_z + rm2_z) / marker_count)];
        
        RDS = vertcat(RDS,rRDS);
    end
   
    if(length(R5M) > length(RDS))
        RDS = vertcat(RDS,zeros(length(R5M)-length(RDS),3));
    else
        RDS = RDS(1:length(R5M),:);
    end
end

%RME
if(reconstruct_arr(2)) 
    %Holds original RME for later graphing
    org_RME = RME;
    RME = [];
    
    [sRTR_x, sRTR_y, sRTR_z] = extract_XYZ(sRTR);
    [sRGT_x, sRGT_y, sRGT_z] = extract_XYZ(sRGT);
    [sRLE_x, sRLE_y, sRLE_z] = extract_XYZ(sRLE);
    [sRME_x, sRME_y, sRME_z] = extract_XYZ(sRME);
    
    %Creates static segment lengths
    sRTR_x = mean(sRTR_x); sRTR_y = mean(sRTR_y); sRTR_z = mean(sRTR_z); 
    sRGT_x = mean(sRGT_x); sRGT_y = mean(sRGT_y); sRGT_z = mean(sRGT_z);
    sRLE_x = mean(sRLE_x); sRLE_y = mean(sRLE_y); sRLE_z = mean(sRLE_z);
    sRME_x = mean(sRME_x); sRME_y = mean(sRME_y); sRME_z = mean(sRME_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(RTR(gc(i):gc(i+1)),RGT(gc(i):gc(i+1)),RLE(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            RME = vertcat(RME,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        %If marker_good is 0, rm1_x,y,z will be 0
        %marker_count and marker_good should be consistent with one another
        marker_good = check_marker(RTR_x(gc(i):gc(i+1)));
        rm1_x = (RTR_x(gc(i):gc(i+1)) + -(dir)*abs(sRTR_x - sRME_x)) * marker_good;
        rm1_y = (RTR_y(gc(i):gc(i+1)) + -(dir)*abs(sRTR_y - sRME_y)) * marker_good;
        rm1_z = (RTR_z(gc(i):gc(i+1)) - abs(sRTR_z - sRME_z)) * marker_good;
        
        marker_good = check_marker(RGT_x(gc(i):gc(i+1)));
        rm2_x = (RGT_x(gc(i):gc(i+1)) + -(dir)*abs(sRGT_x - sRME_x)) * marker_good;
        rm2_y = (RGT_y(gc(i):gc(i+1)) + -(dir)*abs(sRGT_y - sRME_y)) * marker_good;
        rm2_z = (RGT_z(gc(i):gc(i+1)) - abs(sRGT_z - sRME_z)) * marker_good;
        
        marker_good = check_marker(RLE_x(gc(i):gc(i+1)));
        rm3_x = (RLE_x(gc(i):gc(i+1)) + -(dir)*abs(sRLE_x - sRME_x)) * marker_good;
        rm3_y = (RLE_y(gc(i):gc(i+1)) + -(dir)*abs(sRLE_y - sRME_y)) * marker_good;
        rm3_z = (RLE_z(gc(i):gc(i+1)) - abs(sRLE_z - sRME_z)) * marker_good;
        
        rRME = [((rm1_x + rm2_x + rm3_x) / marker_count) ((rm1_y + rm2_y + rm3_y) / marker_count) ((rm1_z + rm2_z + rm3_z) / marker_count)];
        
        RME = vertcat(RME,rRME);   
    end
    
    if(length(R5M) > length(RME))
        RME = vertcat(RME,zeros(length(R5M)-length(RME),3));
    else
        RME = RME(1:length(R5M),:);
    end
end

%R2M
if(reconstruct_arr(6))    
    org_R2M = R2M;
    R2M = [];
    
    [sACB_x, sACB_y, sACB_z] = extract_XYZ(sACB);
    [sR5M_x, sR5M_y, sR5M_z] = extract_XYZ(sR5M);
    [sDLMC5_x, sDLMC5_y, sDLMC5_z] = extract_XYZ(sDLMC5);
    [sR2M_x, sR2M_y, sR2M_z] = extract_XYZ(sR2M);
    
    %Creates static segment lengths
    sACB_x = mean(sACB_x); sACB_y = mean(sACB_y); sACB_z = mean(sACB_z); 
    sR5M_x = mean(sR5M_x); sR5M_y = mean(sR5M_y); sR5M_z = mean(sR5M_z);
    sDLMC5_x = mean(sDLMC5_x); sDLMC5_y = mean(sDLMC5_y); sDLMC5_z = mean(sDLMC5_z);
    sR2M_x = mean(sR2M_x); sR2M_y = mean(sR2M_y); sR2M_z = mean(sR2M_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(ACB(gc(i):gc(i+1)),R5M(gc(i):gc(i+1)),DLMC5(gc(i):gc(i+1)))
        
        if(marker_count == 0)
            R2M = vertcat(R2M,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        %If marker_good is 0, rm1_x,y,z will be 0
        %marker_count and marker_good should be consistent with one another
        marker_good = check_marker(ACB_x(gc(i):gc(i+1)));
        rm1_x = (ACB_x(gc(i):gc(i+1)) + -(dir)*abs(sACB_x - sR2M_x)) * marker_good;
        rm1_y = (ACB_y(gc(i):gc(i+1)) + -(dir)*abs(sACB_y - sR2M_y)) * marker_good;
        rm1_z = (ACB_z(gc(i):gc(i+1)) - abs(sACB_z - sR2M_z)) * marker_good;
        
        marker_good = check_marker(R5M_x(gc(i):gc(i+1)));
        rm2_x = (R5M_x(gc(i):gc(i+1)) + -(dir)*abs(sR5M_x - sR2M_x)) * marker_good;
        rm2_y = (R5M_y(gc(i):gc(i+1)) + -(dir)*abs(sR5M_y - sR2M_y)) * marker_good;
        rm2_z = (R5M_z(gc(i):gc(i+1)) - abs(sR5M_z - sR2M_z)) * marker_good;
        
        marker_good = check_marker(DLMC5_x(gc(i):gc(i+1)));
        rm3_x = (DLMC5_x(gc(i):gc(i+1)) + -(dir)*abs(sDLMC5_x - sR2M_x)) * marker_good;
        rm3_y = (DLMC5_y(gc(i):gc(i+1)) + -(dir)*abs(sDLMC5_y - sR2M_y)) * marker_good;
        rm3_z = (DLMC5_z(gc(i):gc(i+1)) - abs(sDLMC5_z - sR2M_z)) * marker_good;
        
        rR2M = [((rm1_x + rm2_x + rm3_x) / marker_count) ((rm1_y + rm2_y + rm3_y) / marker_count) ((rm1_z + rm2_z + rm3_z) / marker_count)];
        
        R2M = vertcat(R2M,rR2M);   
    end
    
    if(length(R5M) > length(R2M))
        R2M = vertcat(R2M,zeros(length(R5M)-length(R2M),3));
    else
        R2M = R2M(1:length(R5M),:);
    end
end

%RMS
if(reconstruct_arr(4))
    org_RMS = RMS;
    RMS = [];
    [sRLS_x, sRLS_y, sRLS_z] = extract_XYZ(sRLS);
    [sRLO_x, sRLO_y, sRLO_z] = extract_XYZ(sRLO);
    [sRCR_x, sRCR_y, sRCR_z] = extract_XYZ(sRCR);
    [sRMS_x, sRMS_y, sRMS_z] = extract_XYZ(sRMS);
    
    sRLS_x = mean(sRLS_x); sRLS_y = mean(sRLS_y); sRLS_z = mean(sRLS_z); 
    sRLO_x = mean(sRLO_x); sRLO_y = mean(sRLO_y); sRLO_z = mean(sRLO_z);
    sRCR_x = mean(sRCR_x); sRCR_y = mean(sRCR_y); sRCR_z = mean(sRCR_z);
    sRMS_x = mean(sRMS_x); sRMS_y = mean(sRMS_y); sRMS_z = mean(sRMS_z);
    
    for i = 1:gait_cycle_count
        marker_count = count_good_markers(RLS(gc(i):gc(i+1)),RLO(gc(i):gc(i+1)),RCR(gc(i):gc(i+1)));
        
        if(marker_count == 0)
            RMS = vertcat(RMS,zeros(gc(i+1)-gc(i),3));
            continue;
        end
        
        marker_good = check_marker(RLS_x(gc(i):gc(i+1)));
        rm1_x = (RLS_x(gc(i):gc(i+1))+ -(dir)*abs(sRLS_x - sRMS_x)) * marker_good;
        rm1_y = (RLS_y(gc(i):gc(i+1)) + -(dir)*abs(sRLS_y - sRMS_y)) * marker_good;
        rm1_z = (RLS_z(gc(i):gc(i+1)) + abs(sRLS_z - sRMS_z)) * marker_good;
        
        marker_good = check_marker(RLO_x(gc(i):gc(i+1)));
        rm2_x = (RLO_x(gc(i):gc(i+1)) + -(dir)*abs(sRLO_x - sRMS_x)) * marker_good;
        rm2_y = (RLO_y(gc(i):gc(i+1)) + -(dir)*abs(sRLO_y - sRMS_y)) * marker_good;
        rm2_z = (RLO_z(gc(i):gc(i+1)) + abs(sRLO_z - sRMS_z)) * marker_good;
        
        marker_good = check_marker(RCR_x(gc(i):gc(i+1)));
        rm3_x = (RCR_x(gc(i):gc(i+1)) + -(dir)*abs(sRCR_x - sRMS_x)) * marker_good;
        rm3_y = (RCR_y(gc(i):gc(i+1)) + -(dir)*abs(sRCR_y - sRMS_y)) * marker_good;
        rm3_z = (RCR_z(gc(i):gc(i+1)) + abs(sRCR_z - sRMS_z)) * marker_good;
        
        rRMS = [((rm1_x + rm2_x + rm3_x) / marker_count) ((rm1_y + rm2_y + rm3_y) / marker_count) ((rm1_z + rm2_z + rm3_z) / marker_count)];
        
        RMS = vertcat(RMS,rRMS);
    end
   
    if(length(R5M) > length(RMS))
        RMS = vertcat(RMS,zeros(length(R5M)-length(RMS),3));
    else
        RMS = RMS(1:length(R5M),:);
    end
end

if(show_reconst_graphs)   
    %RME
    if(reconstruct_arr(2))
        figure(1)
        plot(new_new_time,org_RME(:,1))
        hold on;
        plot(new_new_time,RTR(:,1))
        plot(new_new_time,RLE(:,1))
        plot(new_new_time,RGT(:,1))
        plot(new_new_time,RME(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed RME data')
        legend('original RME','RTR','RLE','RGT','Reconstructed','Location','Northeast')

        figure(2)
        plot(new_new_time,org_RME(:,2))
        hold on;
        plot(new_new_time,RTR(:,2))
        plot(new_new_time,RLE(:,2))
        plot(new_new_time,RGT(:,2))
        plot(new_new_time,RME(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed RME data')
        legend('original RME','RTR','RLE','RGT','Reconstructed','Location','Northeast')

        figure(3)
        plot(new_new_time,org_RME(:,3))
        hold on;
        plot(new_new_time,RTR(:,3))
        plot(new_new_time,RLE(:,3))
        plot(new_new_time,RGT(:,3))
        plot(new_new_time,RME(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed RME data')
        legend('original RME','RTR','RLE','RGT','Reconstructed','Location','Northeast')
    end
    
    %DS
    if(reconstruct_arr(3))
        figure(19)
        plot(new_new_time,org_RDS(:,1))
        hold on;
        plot(new_new_time,Centroid(:,1))
        plot(new_new_time,RDS(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed RDS data')
        legend('original RDS','Centroid''Reconstructed','Location','Northeast')

        figure(20)
        plot(new_new_time,org_RDS(:,2))
        hold on;
        plot(new_new_time,Centroid(:,2))
        plot(new_new_time,RDS(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed RDS data')
        legend('original RDS','Centroid','Reconstructed','Location','Northeast')

        figure(21)
        plot(new_new_time,org_RDS(:,3))
        hold on;
        plot(new_new_time,Centroid(:,3))
        plot(new_new_time,RDS(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed RDS data')
        legend('original RDS','Centroid','Reconstructed','Location','Northeast')
    end
    
    %RMS
    if(reconstruct_arr(4))
        figure(13)
        plot(new_new_time,org_RMS(:,1))
        hold on;
        plot(new_new_time,RLS(:,1))
        plot(new_new_time,RLO(:,1))
        plot(new_new_time,RCR(:,1))
        plot(new_new_time,RMS(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed RMS data')
        legend('original RMS','RLS','RLO','RCR','Reconstructed','Location','Northeast')

        figure(14)
        plot(new_new_time,org_RMS(:,2))
        hold on;
        plot(new_new_time,RLS(:,2))
        plot(new_new_time,RLO(:,2))
        plot(new_new_time,RCR(:,2))
        plot(new_new_time,RMS(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed RMS data')
        legend('original RMS','RLS','RLO','RCR','Reconstructed','Location','Northeast')

        figure(15)
        plot(new_new_time,org_RMS(:,3))
        hold on;
        plot(new_new_time,RLS(:,3))
        plot(new_new_time,RLO(:,3))
        plot(new_new_time,RCR(:,3))
        plot(new_new_time,RMS(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed RMS data')
        legend('original RMS','RLS','RLO','RCR','Reconstructed','Location','Northeast')
    end
    
   
    
    if(reconstruct_arr(5))
        figure(4)
        plot(new_new_time,org_ACB(:,1))
        hold on;
        plot(new_new_time,R5M(:,1))
        plot(new_new_time,R2M(:,1))
        plot(new_new_time,DLMC5(:,1))
        plot(new_new_time,ACB(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed ACB data')
        legend('original ACB','R5M','R2M','DLMC5','Reconstructed','Location','Northeast')

        figure(5)
        plot(new_new_time,org_ACB(:,2))
        hold on;
        plot(new_new_time,R5M(:,2))
        plot(new_new_time,R2M(:,2))
        plot(new_new_time,DLMC5(:,2))
        plot(new_new_time,ACB(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed ACB data')
        legend('original ACB','R5M','R2M','DLMC5','Reconstructed','Location','Northeast')

        figure(6)
        plot(new_new_time,org_ACB(:,3))
        hold on;
        plot(new_new_time,R5M(:,3))
        plot(new_new_time,R2M(:,3))
        plot(new_new_time,DLMC5(:,3))
        plot(new_new_time,ACB(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed ACB data')
        legend('original ACB','R5M','R2M','DLMC5','Reconstructed','Location','Northeast')
    end
    
    %MP2
    if(reconstruct_arr(6))
        figure(22)
        plot(new_new_time,org_R2M(:,1))
        hold on;
        plot(new_new_time,ACB(:,1))
        plot(new_new_time,R5M(:,1))
        plot(new_new_time,DLMC5(:,1))
        plot(new_new_time,R2M(:,1))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('x values of reconstructed R2M data')
        legend('original R2M','ACB','R5M','DLMC5','Reconstructed','Location','Northeast')

        figure(23)
        plot(new_new_time,org_R2M(:,2))
        hold on;
        plot(new_new_time,ACB(:,2))
        plot(new_new_time,R5M(:,2))
        plot(new_new_time,DLMC5(:,2))
        plot(new_new_time,R2M(:,2))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('y values of reconstructed R2M data')
        legend('original R2M','ACB','R5M','DLMC5','Reconstructed','Location','Northeast')

        figure(24)
        plot(new_new_time,org_R2M(:,3))
        hold on;
        plot(new_new_time,ACB(:,3))
        plot(new_new_time,R5M(:,3))
        plot(new_new_time,DLMC5(:,3))
        plot(new_new_time,R2M(:,3))
        hold off;
        xlabel('time(s)')
        ylabel('pos (mm)')
        title('z values of reconstructed R2M data')
        legend('original R2M','ACB','R5M','DLMC5','Reconstructed','Location','Northeast')
    end
 
end

oc = 1;

%Error checking
for n = 1:gait_cycle_count
   
    start_frame = gc(n);
    end_frame = gc(n+1);
    segment_error_checks = struct("RGT_RLE",0,"RLO_RLS",0,"RSC1_RAC",0,"ACB_R5M",0);
    
    segment_error_checks.RGT_RLE = error_check(RGT,RLE,error_length,static_RGT_RLE,start_frame,end_frame,1);
    
    if(override_arr(oc))
        segment_error_checks.RGT_RLE = ~segment_error_checks.RGT_RLE;
    end
    
    segment_error_checks.RLO_RLS = error_check(RLO,RLS,error_length,static_RLO_RLS,start_frame,end_frame,0);
    
    if(override_arr(oc + 1))
        segment_error_checks.RLO_RLS = ~segment_error_checks.RLO_RLS;
    end
    
    segment_error_checks.RSC1_RAC = error_check(RSC1,RAC,error_length,static_RSC1_RAC,start_frame,end_frame,0);
    
    if(override_arr(oc + 2))
        segment_error_checks.RSC1_RAC = ~segment_error_checks.RSC1_RAC;
    end
    
    segment_error_checks.ACB_R5M = error_check(ACB,R5M,error_length,static_ACB_R5M,start_frame,end_frame);

    if(override_arr(oc + 3))  
        segment_error_checks.ACB_R5M = ~segment_error_checks.ACB_R5M;
    end
          
    gait_cycles(n) = struct("start_frame",start_frame,"end_frame",end_frame,"seg_checks",segment_error_checks);
    
    oc = oc + 4;
end

if(show_graphs) 
    rgt_rle = lengthPlotter(static_RGT_RLE,RGT,RLE);
    new_new_time = 1:length(rgt_rle);
    
    figure(8)
    plot(new_new_time,rgt_rle)
    hold on;
    yline(static_RGT_RLE);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Humerus Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')
   
    rlo_rls = lengthPlotter(static_RLO_RLS,RLO,RLS);
    figure(9)
    plot(new_new_time,rlo_rls)
    hold on;
    yline(static_RLO_RLS);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Radius/Ulna Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')

    rds_cent = lengthPlotter(static_RSC1_RAC,RDS,Centroid);
    figure(10)
    plot(new_new_time,rds_cent)
    hold on;
    yline(static_RSC1_RAC);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Scalpula Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')

    acb_r5m = lengthPlotter(static_ACB_R5M,ACB,R5M);
    figure(12)
    plot(new_new_time,acb_r5m)
    hold on;
    yline(static_ACB_R5M);
    hold off;
    xlabel('frames')
    ylabel('length(mm)')
    title('Length of Static and Dynamic Manus Segments vs Time')
    legend('Dynamic segment','Static segment','Location','Northwest')
end

trial = struct("trial_name",-1,"gait_cycles",-1);

name = extractBetween(name,1,length(name)-4);
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
    tACB = ACB(start_frame:end_frame,:);
    
    x(i) = struct("is_good",gc_is_good,"R5M", tR5M, "R2M", tR2M, "RGT", tRGT, "RLE", tRLE, ...
    "RLO", tRLO, "RLS", tRLS, "T1", tT1, "RDS", tRDS, "Centroid", tCentroid, "RME", tRME, ...
    "RTR", tRTR, "RMS", tRMS, "DLMC5", tDLMC5, "ACB", tACB, "start_frame", start_frame, "end_frame",end_frame, ...
    "RGT_RLE", gait_cycles(i).seg_checks.RGT_RLE, "RLO_RLS", gait_cycles(i).seg_checks.RLO_RLS, ...
    "RSC1_RAC", gait_cycles(i).seg_checks.RSC1_RAC, ...
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
trial.pos_data.RSC1 = RSC1;
trial.pos_data.RCDS = RCDS;
trial.pos_data.RSC1 = RSC1;
trial.pos_data.RAC = RAC;

trial.gait_cycles = x;

save(['Produced Data/' name '.mat'],'trial')
toc


end
      