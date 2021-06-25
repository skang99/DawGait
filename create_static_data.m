function [static_RGT_RLE,static_RLO_RLS,static_RDS_Centroid,static_T1_Centroid,static_ACB_R5M,T1,VTR1,SC1,DLMC5,ACB,RDS,RTR,RME,RCR,R5M,R2M,RMS,RGT,RLE,Centroid,RLS,RLO] = create_static_data(static_trial)

try 
[static_data,static_coords] = xlsread(static_trial); %import excel data numeric and text
catch exception
    disp("Error in opening the specified static Excel file. File may not exist or is not within the current directory.");
    return;
end

%New system
static_pos_data = static_data(6:length(static_data),:); %Extract only the numeric data a.k.a position readings

RGT = double(subs(extract_data(static_pos_data,static_coords,'GRTB'),NaN,0));
RLE = double(subs(extract_data(static_pos_data,static_coords,'LEPI'),NaN,0));
R5M = double(subs(extract_data(static_pos_data,static_coords,'MCP5'),NaN,0));
R2M = double(subs(extract_data(static_pos_data,static_coords,'MCP2'),NaN,0));
RLO = double(subs(extract_data(static_pos_data,static_coords,'OLEC'),NaN,0));
RLS = double(subs(extract_data(static_pos_data,static_coords,'LSTY'),NaN,0));
T1 = double(subs(extract_data(static_pos_data,static_coords,'T1'),NaN,0));
RDS = double(subs(extract_data(static_pos_data,static_coords,'DRSC'),NaN,0));
RME = double(subs(extract_data(static_pos_data,static_coords,'MEPI'),NaN,0));
RTR = double(subs(extract_data(static_pos_data,static_coords,'TRIC'),NaN,0));
RCR = double(subs(extract_data(static_pos_data,static_coords,'MDCR'),NaN,0));
RMS = double(subs(extract_data(static_pos_data,static_coords,'MSTY'),NaN,0));

RAC = double(subs(extract_data(static_pos_data,static_coords,'ACRM'),NaN,0));
RSC1 = double(subs(extract_data(static_pos_data,static_coords,'SC1'),NaN,0));
RSC2 = double(subs(extract_data(static_pos_data,static_coords,'SC2'),NaN,0));

VTR1 = double(subs(extract_data(static_pos_data,static_coords,'VTR1'),NaN,0));
SC1 = double(subs(extract_data(static_pos_data,static_coords,'SC1'),NaN,0));
DLMC5 = double(subs(extract_data(static_pos_data,static_coords,'DLMC5'),NaN,0));
ACB = double(subs(extract_data(static_pos_data,static_coords,'ACCB'),NaN,0));

% static_pos_data = static_data(7:length(static_data),:); %Extract only the numeric data a.k.a position readings
%  
% RGT = extract_data(static_pos_data,static_coords,'R Greater Tubercle');
% RLE = extract_data(static_pos_data,static_coords,'R Lateral Epicondyle');
% R5M = extract_data(static_pos_data,static_coords,'R 5th Metacarpal');
% RLO = extract_data(static_pos_data,static_coords,'R Lateral Olecranon');
% RLS = extract_data(static_pos_data,static_coords,'R Lateral Styloid');
% T1 = extract_data(static_pos_data,static_coords,'T1');
% RDS = extract_data(static_pos_data,static_coords,'R Dorsal Scapula');
%  
% RAC = extract_data(static_pos_data,static_coords,'R Acromion');
% RSC1 = extract_data(static_pos_data,static_coords,'RSC1');
% RSC2 = extract_data(static_pos_data,static_coords,'RSC2');
% 
% RME = extract_data(static_pos_data,static_coords,'R Medial Epicondyle');
% RTR = extract_data(static_pos_data,static_coords,'R Tricep');
% RMS = extract_data(static_pos_data,static_coords,'R Medial Styloid');
 

Centroid = (RAC + RSC1 + RSC2) / 3;

R_5th_M_z = RGT(:,3);

[t_naught, t_one] = find_still_static(R_5th_M_z);

%Current work around. Frame # starts at 15, 46, etc, but the start should
%be t = 0, which this does not currently account for

% lower_lim =  %finds the frame where the lower limit time value is located
% upper_lim = find(time == t_one) %finds the frame where the upper limit time value is located

lower_lim =  1; %finds the frame where the lower limit time value is located
upper_lim =  100;%finds the frame where the upper limit time value is located

RGT = RGT(lower_lim: upper_lim,:); %extracts each frame of data within the specified range
RLE = RLE(lower_lim: upper_lim,:); 
R5M = R5M(lower_lim: upper_lim,:);
RLO = RLO(lower_lim: upper_lim,:); 
RLS = RLS(lower_lim: upper_lim,:); 
T1 = T1(lower_lim: upper_lim,:);
RDS = RDS(lower_lim: upper_lim,:);
Centroid = Centroid(lower_lim: upper_lim,:);
ACB = ACB(lower_lim:upper_lim,:);
% RME = RME(lower_lim: upper_lim,:);
% RTR = RTR(lower_lim: upper_lim,:);


ab_x = RGT(:,1) - RLE(:,1); 
ab_z = RGT(:,3) - RLE(:,3); 
ab = sqrt(ab_x.^2 + ab_z.^2);

%Averages the previous distances to create a scalar. 
%This scalar value will be used as a baseline to
%determine whether or not dynamic gait cycles
%have too much skin movement in the error_check function.
static_RGT_RLE = mean(ab); 
                      
de_x = RLO(:,1) - RLS(:,1);
de_z = RLO(:,3) - RLS(:,3);
de = sqrt(de_x.^2 + de_z.^2); 
static_RLO_RLS = mean(de); 

gj_x = RDS(:,1) - Centroid(:,1); 
gj_z = RDS(:,3) - Centroid(:,3); 
gj = sqrt(gj_x.^2 + gj_z.^2); 
static_RDS_Centroid = mean(gj);                       
                           
hb_x = T1(:,1) - Centroid(:,1); 
hb_z = T1(:,3) - Centroid(:,3); 
hb = sqrt(hb_x.^2 + hb_z.^2);
static_T1_Centroid = mean(hb);

am_x = ACB(:,1) - R5M(:,1);
am_z = ACB(:,3) - R5M(:,3);
am = sqrt(am_x.^2 + am_z.^2);
static_ACB_R5M = mean(am);



                                                                                
end







